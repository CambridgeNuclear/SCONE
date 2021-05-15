module asciiOutput_inter

  use numPrecision

  implicit none
  private

  !!
  !! Abstract interface for output printers
  !!
  !! Allows stream-like output to multiple formats (e.g. MATLAB, CSV, JSON)
  !!
  !! They recive the sequence of calls from the `outputFile` and conver them to the output file.
  !!
  !! Printers do no error checking. It is responsibility of the `outputFile` to ensure that
  !! the sequence is correct.
  !!
  !! NOTE:
  !!  outputFile converts the numerical values to characters. Thus printer gets reals already as
  !!  chars.
  !!
  !! Interface:
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
  contains
    procedure(writeToFile),deferred :: writeToFile
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
    !! Print the output to the given unit
    !!
    !! Args:
    !!  unit [in] -> Unit number of the output
    !!
    subroutine writeToFile(self, unit)
      import :: asciiOutput, &
                shortInt
      class(asciiOutput), intent(inout) :: self
      integer(shortInt), intent(in)     :: unit
    end subroutine writeToFile

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
    !! Can recive single value or array next
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
    !! Name should alrady be provided by "startEntry"
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

end module asciiOutput_inter
