module asciiJSON_class

  use numPrecision
  use genericProcedures, only : numToChar
  use asciiOutput_inter, only : asciiOutput
  use charTape_class,    only : charTape

  implicit none
  private

  character(1), parameter :: BLANK = ' ', &
                             QUOTE = '"', &
                             NEWLINE = new_line(BLANK)
  integer(shortInt), parameter :: IND_SIZE = 2

  !!
  !! Printer for JSON output
  !!
  !! Creates a JSON file with an output
  !! It is in ASCII and is semi-human readable.
  !! There is indentation for each block and entry, but long multi-dimensional (N-D) arrays
  !! are printed on a single line. Also N-D arrays are represented as a JSON array of arrays.
  !!
  !! JSON output is fairly standard so it can be easily read into Python and other environments
  !!
  !! TODO:
  !!  Currently there is a dirty fix to remove commas at the end of a block, which requires to
  !!  rewind the output. It is not ideal and could be improved.
  !!
  !! Private Members:
  !!   output      -> Character that contain the output
  !!   ind_lvl     -> Indentation level counter
  !!   shapeBuffer -> Shape of the current array
  !!   in_array    -> True is is writing an array
  !!   count       -> Number of entries written to the array
  !!
  !! Interface:
  !!   asciiOutput Interface
  !!
  type, public, extends(asciiOutput) :: asciiJSON
    private
    integer(shortInt) :: ind_lvl = 0

    integer(shortInt), dimension(:), allocatable :: shapeBuffer
    logical(defBool)                             :: in_array = .false.
    integer(shortInt)                            :: count = 0


  contains
    procedure :: init
    procedure :: endFile
    procedure :: extension

    procedure :: startBlock
    procedure :: endBlock
    procedure :: startEntry
    procedure :: endEntry
    procedure :: startArray
    procedure :: endArray
    procedure :: printNum
    procedure :: printChar
  end type asciiJSON

contains

  !!
  !! Initialise the JSON output
  !!
  subroutine init(self)
    class(asciiJSON), intent(inout) :: self

    ! Add the initial bracket
    call self % append( "{" // NEWLINE)
    self % ind_lvl = 1

  end subroutine init

  !!
  !! Finalise the printer
  !!
  !! See asciiOutput_inter for details
  !!
  subroutine endFile(self)
    class(asciiJSON), intent(inout) :: self

    ! Dirty fix
    ! Remove comma from the last entry by rewind
    ! Need to check that the character is a comma to avoid removing { when there is an empty block
    ! TODO: Find better way. Will clash with any proper stream output
    if (self % peek(2) == ",") then
      call self % cut(2)
      call self % append(NEWLINE)
    end if
    call self % append("}" // NEWLINE)

  end subroutine endFile

  !!
  !! Return appropriate extension for the file
  !!
  !! See asciiOutput_inter for details
  !!
  pure function extension(self) result(str)
    class(asciiJSON), intent(in) :: self
    character(:), allocatable    :: str

    str = 'json'

  end function extension

  !!
  !! Change state to writing new block with "name"
  !!
  !! See asciiOutput_inter for details
  !!
  subroutine startBlock(self, name)
    class(asciiJSON), intent(inout) :: self
    character(nameLen), intent(in)  :: name

    ! Write indentation
    call self % append(repeat(BLANK, self % ind_lvl * IND_SIZE))

    ! Open new block and increase indentation
    call self % append(QUOTE // trim(name) // QUOTE // ":" // "{" // NEWLINE)
    self % ind_lvl = self % ind_lvl + 1

  end subroutine startBlock

  !!
  !! End current block
  !!
  !! See asciiOutput_inter for details
  !!
  subroutine endBlock(self)
    class(asciiJSON), intent(inout) :: self

    ! Dirty fix
    ! Remove comma from the last entry by rewind
    ! Need to check that the character is a comma to avoid removing { when there is an empty block
    ! TODO: Find better way. Will clash with any proper stream output
    if (self % peek(2)  == ",") then
      call self % cut(2)
      call self % append(NEWLINE)
    end if

    ! Decrease and write indentation
    self % ind_lvl = self % ind_lvl - 1
    call self % append(repeat(BLANK, self % ind_lvl * IND_SIZE))

    ! Close the block
    call self % append("}" // ","// NEWLINE)

  end subroutine endBlock

  !!
  !! Change state to writing a new entry
  !!
  !! See asciiOutput_inter for details
  !!
  subroutine startEntry(self, name)
    class(asciiJSON), intent(inout) :: self
    character(*), intent(in)        :: name

    ! Print indentation
    call self % append(repeat(BLANK, self % ind_lvl * IND_SIZE))

    ! Write the name
    call self % append(QUOTE // trim(name) // QUOTE // ":")

  end subroutine startEntry

  !!
  !! End writing a new entry
  !!
  !! See asciiOutput_inter for details
  !!
  subroutine endEntry(self)
    class(asciiJSON), intent(inout) :: self

    ! End with NEWLINE
    ! Comma will be written together with a number/char
    call self % append(NEWLINE)

  end subroutine endEntry

  !!
  !! Start writing array with the given shape
  !!
  !! See asciiOutput_inter for details
  !!
  subroutine startArray(self, shape)
    class(asciiJSON), intent(inout)           :: self
    integer(shortInt),dimension(:),intent(in) :: shape

    ! Save the shape buffer
    self % shapeBuffer = shape
    self % in_array = .true.
    self % count = 0

  end subroutine startArray

  !!
  !! End writing the array
  !!
  !! See asciiOutput_inter for details
  !!
  subroutine endArray(self)
    class(asciiJSON), intent(inout) :: self

    self % in_array = .false.

  end subroutine endArray

  !!
  !! Print numerical value
  !!
  !! See asciiOutput_inter for details
  !!
  subroutine printNum(self, val)
    class(asciiJSON), intent(inout) :: self
    character(*),intent(in)         :: val
    integer(shortInt)               :: i, mul

    if (self % in_array) then
      ! Print opening brackets if required
      mul = 1
      do i = 1, size(self % shapeBuffer)
        mul = mul * self % shapeBuffer(i)
        if (modulo(self % count, mul) == 0) call self % append("[")
      end do
    end if

    ! Print the number
    call self % append(val)

    if (self % in_array) then
      ! Append the count
      self % count = self % count + 1

      ! Print required brackets
      mul = 1
      do i = 1, size(self % shapeBuffer)
        mul = mul * self % shapeBuffer(i)
        if (modulo(self % count, mul) == 0) call self % append("]")
      end do
    end if

    ! Finish entry with a comma
    call self % append(",")

  end subroutine printNum

  !!
  !! Print character value
  !!
  !! See asciiOutput_inter for details
  !!
  subroutine printChar(self, val)
    class(asciiJSON), intent(inout) :: self
    character(*),intent(in)         :: val
    integer(shortInt)               :: i, mul

    if (self % in_array) then
      ! Print opening brackets if required
      mul = 1
      do i = 1, size(self % shapeBuffer)
        mul = mul * self % shapeBuffer(i)
        if (modulo(self % count, mul) == 0) call self % append("[")
      end do
    end if

    ! Print the character
    call self % append(QUOTE // val // QUOTE)

    if (self % in_array) then
      ! Append the count
      self % count = self % count + 1

      ! Print required brackets
      mul = 1
      do i = 1, size(self % shapeBuffer)
        mul = mul * self % shapeBuffer(i)
        if (modulo(self % count, mul) == 0) call self % append("]")
      end do
    end if

    ! Finish entry with a comma
    call self % append(",")

  end subroutine printChar

end module asciiJSON_class
