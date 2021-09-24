module asciiMATLAB_class

  use numPrecision
  use genericProcedures, only : numToChar
  use stack_class,       only : stackChar
  use asciiOutput_inter, only : asciiOutput
  use charTape_class,    only : charTape

  implicit none
  private

  !! Local parameters
  character(1),parameter :: BLANK = ' ', &
                            NEWLINE = new_line(BLANK) ,&
                            APOS    = char(39) ,&
                            BRAKET_L  = '{' ,&
                            BRAKET_R  = '}'

  integer(shortInt), parameter :: IN_BLOCK = 0, &
                                  IN_ENTRY = 1, &
                                  IN_ARRAY = 2

  !!
  !! Printer for ASCII MATLAB output file
  !!
  !! Creates a computer-readable matalab ".m" file.
  !! Despite beeing in ASCII it is not intended for human-readability:
  !!  -> Each entry is a single line.
  !!  -> Furthermore every array is printed as RANK 1 array and a reshape function.
  !!
  !! However the above scheme ensures quick reading by MATLAB [no reallocation of arrays.
  !! Looking at you Serpent :P ]
  !!
  !! Each block is marked with a prefix in the variable name e.g.:
  !!   BLOCK1_BLOCK2_Entry1 = reshape([1,2,3,4],2,2);
  !!
  !! Private members:
  !!   state          -> Current state
  !!   blockNameStack -> Names of all blocks to construct the prefix
  !!   blockLevel     -> Current block level
  !!   output         -> Character that contain the output
  !!   prefix         -> Current prefix for the block
  !!   shapeBuffer    -> Buffored shape of the current array
  !!
  !! Interface:
  !!   asciiOutput Interface
  !!
  type, public, extends(asciiOutput) :: asciiMATLAB
    private
    ! State components
    integer(shortInt) :: state  = IN_BLOCK
    type(stackChar)   :: blockNameStack
    integer(shortInt) :: blockLevel = 0

    ! Outputfile
    type(charTape) :: output

    ! Buffers
    type(charTape) :: prefix
    integer(shortInt), dimension(:), allocatable :: shapeBuffer

  contains
    procedure :: init
    procedure :: extension
    procedure :: writeToFile

    procedure :: startBlock
    procedure :: endBlock
    procedure :: startEntry
    procedure :: endEntry
    procedure :: startArray
    procedure :: endArray
    procedure :: printNum
    procedure :: printChar

  end type asciiMATLAB

contains

  !!
  !! Initialise the printer
  !!
  !! See asciiOutput_inter for details
  !!
  subroutine init(self)
    class(asciiMATLAB), intent(inout) :: self

    ! Nothing to do

  end subroutine init

  !!
  !! Return approperiate extension for the file
  !!
  !! See asciiOutput_inter for details
  !!
  pure function extension(self) result(str)
    class(asciiMATLAB), intent(in) :: self
    character(:), allocatable      :: str

    str = 'm'

  end function extension

  !!
  !! Print the output to the given unit
  !!
  !! See asciiOutput_inter for details
  !!
  subroutine writeToFile(self, unit)
    class(asciiMATLAB), intent(inout) :: self
    integer(shortInt), intent(in)     :: unit
    character(:),allocatable          :: form

    form = '(A' // numToChar(self % output % length()) // ')'

    write(unit,form) self % output % expose()

  end subroutine writeToFile

  !!
  !! Change state to writing new block with "name"
  !!
  !! See asciiOutput_inter for details
  !!
  subroutine startBlock(self, name)
    class(asciiMATLAB), intent(inout) :: self
    character(nameLen), intent(in)    :: name
    character(100), parameter :: Here ='startBlock (asciiMATLAB_class.f90)'

    ! Check if state support acction
    ! if (self % state == IN_ENTRY .or. self % state == IN_ARRAY) then
    !   call fatalError(Here,'Cannot start writing new block inside entry or array')
    ! end if

    ! Update state - change current prefix
    call self % blockNameStack % push(name)
    call self % prefix % append(trim(name) // '_')
    self % blockLevel = self % blockLevel + 1

  end subroutine startBlock

  !!
  !! End current block
  !!
  !! See asciiOutput_inter for details
  !!
  subroutine endBlock(self)
    class(asciiMATLAB), intent(inout) :: self
    character(nameLen)                :: temp
    integer(shortInt)                 :: N
    character(100), parameter :: Here ='endBlock (asciiMATLAB_class.f90)'

    ! Check if state support acction
    ! if (self % state == IN_ENTRY .or. self % state == IN_ARRAY) then
    !   call fatalError(Here,'Cannot end writing new block inside entry or array')
    ! end if
    !
    ! if ( self % blockLevel == 0) then
    !   call fatalError(Here,'Cannot exit from root block')
    ! end if

    ! Update state - change current prefix
    call self % blockNameStack % pop(temp)

    ! Calculate how much prefix needs to be cut
    N = len_trim(temp) + 1
    call self % prefix % cut(N)
    self % blockLevel = self % blockLevel - 1

  end subroutine endBlock

  !!
  !! Change state to writing a new entry
  !!
  !! See asciiOutput_inter for details
  !!
  subroutine startEntry(self, name)
    class(asciiMATLAB), intent(inout) :: self
    character(*), intent(in)          :: name
    character(100), parameter :: Here ='startEntry (asciiMATLAB_class.f90)'

    ! Check if state support acction
    ! if (self % state == IN_ENTRY .or. self % state == IN_ARRAY) then
    !   call fatalError(Here,'Cannot star writing new entry inside entry or array')
    ! end if

    ! Update state
    self % state = IN_ENTRY

    ! Write variable name with prefix
    call self % output % append( trim(self % prefix % expose()) // trim(name))
    call self % output % append( ' = ')


  end subroutine startEntry

  !!
  !! End writing a new entry
  !!
  !! See asciiOutput_inter for details
  !!
  subroutine endEntry(self)
    class(asciiMATLAB), intent(inout) :: self
    character(100), parameter :: Here ='endEntry (asciiMATLAB_class.f90)'

    ! Check if state support acction
    ! if (self % state == IN_BLOCK .or. self % state == IN_ARRAY) then
    !   call fatalError(Here,'Cannot finish entry in block or array')
    ! end if

    ! Update state
    self % state = IN_BLOCK

    ! Write semicolon and newline
    call self % output % append( ';' // NEWLINE)

  end subroutine endEntry

  !!
  !! Start writing array with the given shape
  !!
  !! See asciiOutput_inter for details
  !!
  subroutine startArray(self, shape)
    class(asciiMATLAB), intent(inout)         :: self
    integer(shortInt),dimension(:),intent(in) :: shape
    character(100), parameter :: Here ='startArray (asciiMATLAB_class.f90)'

    ! Check if state support acction
    ! if (self % state /= IN_ENTRY) then
    !   call fatalError(Here,'Cannot finish entry in block or array')
    ! end if

    ! Update state
    self % state = IN_ARRAY

    ! Store shape information
    self % shapeBuffer = shape

    ! Write start of array
    if( size(self % shapeBuffer) == 1) then
      call self % output % append('[ ')

    else
      call self % output % append('reshape([ ')

    end if

  end subroutine startArray

  !!
  !! End writing the array
  !!
  !! See asciiOutput_inter for details
  !!
  subroutine endArray(self)
    class(asciiMATLAB), intent(inout) :: self
    integer(shortInt)                 :: i
    character(100), parameter :: Here ='endArray (asciiMATLAB_class.f90)'

    ! Check if state support acction
    ! if (self % state /= IN_ARRAY) then
    !   call fatalError(Here,'Cannot finish array inside an entry')
    ! end if

    ! Update state
    self % state = IN_ENTRY

    ! Write end of array
    call self % output % cut(1)
    call self % output % append(']')

    ! Finish reshape function for higher rank arrays
    if( size(self % shapeBuffer) > 1) then
      do i=1,size(self % shapeBuffer)
        call self % output % append(',' // numToChar(self % shapeBuffer(i)))
      end do
      call self % output % append(')')
    end if

  end subroutine endArray

  !!
  !! Print numerical value
  !!
  !! See asciiOutput_inter for details
  !!
  subroutine printNum(self, val)
    class(asciiMATLAB), intent(inout) :: self
    character(*),intent(in)           :: val
    character(100), parameter :: Here ='printNum (asciiMATLAB_class.f90)'

    if(self % state == IN_ARRAY) then
      call self % output % append(val//',')

    else if( self % state == IN_ENTRY) then
      call self % output % append(val)

    ! else
    !   call fatalError(Here,'Cannot print number directly into block')
    end if

  end subroutine printNum

  !!
  !! Print character value
  !!
  !! See asciiOutput_inter for details
  !!
  subroutine printChar(self, val)
    class(asciiMATLAB), intent(inout) :: self
    character(*),intent(in)           :: val
    character(100), parameter :: Here ='printChar (asciiMATLAB_class.f90)'

    if(self % state == IN_ARRAY) then
      call self % output % append(BRAKET_L // APOS // val // APOS // BRAKET_R // ",")

    else if( self % state == IN_ENTRY) then
      call self % output % append(BRAKET_L // APOS // val // APOS // BRAKET_R )

    ! else
    !   call fatalError(Here,'Cannot print number directly into block')
    end if

  end subroutine printChar

end module asciiMATLAB_class
