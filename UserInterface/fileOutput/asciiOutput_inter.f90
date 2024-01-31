module asciiOutput_inter

  use numPrecision
  use delayedStream_class, only : delayedStream

  implicit none
  private

  !!
  !! Abstract interface for output printers
  !!
  !! Allows stream-like output to multiple formats (e.g. MATLAB, CSV, JSON)
  !!
  !! They receive the sequence of calls from the `outputFile` and convert them to the output file.
  !!
  !! Printers do no error checking. It is responsibility of the `outputFile` to ensure that
  !! the sequence is correct.
  !!
  !! NOTE:
  !!  outputFile converts the numerical values to characters. Thus printer gets reals already as
  !!  chars.
  !!
  !! Interface:
  !!   init        -> Initialise
  !!   extension   -> Return file extension appropriate for the format
  !!   writeToFile -> Print the output to the provided unit
  !!   startBlock  -> Start new block
  !!   endBlock    -> End a block
  !!   startArray  -> Start array with a given shape
  !!   endArray    -> End current array
  !!   printNum    -> Print a number (given as character)
  !!   printChar   -> Print a character
  !!
  type, public,abstract :: asciiOutput
    private
    type(delayedStream) :: stream
  contains
    procedure                       :: append
    procedure                       :: cut
    procedure                       :: peek
    procedure                       :: close
    procedure                       :: setUnit
    procedure(init), deferred       :: init
    procedure(init), deferred       :: endFile
    procedure(extension), deferred  :: extension
    procedure(startBlock),deferred  :: startBlock
    procedure(endBlock),deferred    :: endBlock
    procedure(startEntry),deferred  :: startEntry
    procedure(endEntry),deferred    :: endEntry
    procedure(startArray),deferred  :: startArray
    procedure(endArray),deferred    :: endArray
    procedure(printNum),deferred    :: printNum
    procedure(printChar),deferred   :: printChar
  end type asciiOutput


  abstract interface

    !!
    !! Initialise the printer
    !!
    !! Args:
    !!  None
    !!
    subroutine init(self)
      import :: asciiOutput
      class(asciiOutput), intent(inout) :: self
    end subroutine init

    !!
    !! End the file
    !!
    !! Format may require to print some characters at the end, so the printer
    !! needs to be notified that no more data will be provided.
    !!
    !! Args:
    !!  None
    !!
    subroutine endFile(self)
      import :: asciiOutput
      class(asciiOutput), intent(inout) :: self
    end subroutine endFile

    !!
    !! Return appropriate extension for the file
    !!
    !! Must be without any "." Thus "exe" instead of ".exe"!
    !!
    !! Args:
    !!   None
    !!
    pure function extension(self) result(str)
      import :: asciiOutput
      class(asciiOutput), intent(in) :: self
      character(:), allocatable      :: str
    end function extension

    !!
    !! Change state to writing new block with "name"
    !!
    !! Args:
    !!   name [in] -> Name of the new block
    !!
    subroutine startBlock(self, name)
      import :: asciiOutput, &
                nameLen
      class(asciiOutput), intent(inout) :: self
      character(nameLen), intent(in)    :: name
    end subroutine startBlock

    !!
    !! End current block
    !!
    subroutine endBlock(self)
      import :: asciiOutput
      class(asciiOutput), intent(inout) :: self
    end subroutine endBlock

    !!
    !! Change state to writing a new entry
    !!
    !! Can receive single value or array next
    !!
    !! Args:
    !!   name [in] -> name of the entry
    !!
    subroutine startEntry(self, name)
      import :: asciiOutput, &
                nameLen
      class(asciiOutput), intent(inout) :: self
      character(nameLen), intent(in)    :: name
    end subroutine startEntry

    !!
    !! End writing a new entry
    !!
    subroutine endEntry(self)
      import :: asciiOutput
      class(asciiOutput), intent(inout) :: self
    end subroutine endEntry

    !!
    !! Start writing array with the given shape
    !!
    !! Name should already be provided by "startEntry"
    !!
    !! Args:
    !!   shape [in] -> Shape of the array (in column-major order)
    !!
    subroutine startArray(self, shape)
      import :: asciiOutput, &
                shortInt
      class(asciiOutput), intent(inout)         :: self
      integer(shortInt),dimension(:),intent(in) :: shape
    end subroutine startArray

    !!
    !! End writing the array
    !!
    subroutine endArray(self)
      import :: asciiOutput
      class(asciiOutput), intent(inout) :: self
    end subroutine endArray

    !!
    !! Print numerical value
    !!
    !! Assume it contains valid number with NO leading or trailing blanks
    !!
    !! Args:
    !!  val [in] -> A number converted to character
    !!
    subroutine printNum(self, val)
      import :: asciiOutput
      class(asciiOutput), intent(inout) :: self
      character(*),intent(in)           :: val
    end subroutine printNum

    !!
    !! Print character value
    !!
    !! Assuming it contains NO leading or trailing blanks
    !!
    !! Args:
    !!  val [in] -> Character value
    !!
    subroutine printChar(self, val)
      import :: asciiOutput
      class(asciiOutput), intent(inout) :: self
      character(*),intent(in)           :: val
    end subroutine printChar
  end interface

contains


    !!
    !! Write data to the output
    !!
    !! Args:
    !!  text [in] -> Text to be written
    !!
    subroutine append(self, text)
      class(asciiOutput), intent(inout) :: self
      character(*), intent(in)          :: text

      call self % stream % write(text)

    end subroutine append

    !!
    !! Remove n characters from the output
    !!
    !! Args:
    !!  n [in] -> Number of characters to be removed
    !!
    subroutine cut(self, n)
      class(asciiOutput), intent(inout) :: self
      integer, intent(in)               :: n

      call self % stream % cut(n)

    end subroutine cut

    !!
    !! Get the nth last character from the output
    !!
    !! Args:
    !!   n [in] -> How many characters back to peak
    !!
    function peek(self, n) result(c)
      class(asciiOutput), intent(inout) :: self
      integer, intent(in)               :: n
      character(1)                      :: c

      c = self % stream % peek(n)

    end function peek

    !!
    !! Close the stream
    !!
    !! Signals to the output that no more data will be provided, so the
    !! closing characters can be printed. In addition flushes remaining characters
    !! in the buffer to the output.
    !!
    subroutine close(self)
      class(asciiOutput), intent(inout) :: self

      call self % endFile()
      call self % stream % close()

    end subroutine close

    !!
    !! Set the unit to write to
    !!
    subroutine setUnit(self, unit)
      class(asciiOutput), intent(inout) :: self
      integer(shortInt), intent(in)     :: unit

      call self % stream % setUnit(unit)

    end subroutine setUnit


end module asciiOutput_inter
