module asciiMATLAB_class

  use numPrecision
  use genericProcedures, only : fatalError, numToChar
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
  !! Constructor
  !!
  interface asciiMATLAB
    module procedure asciiMATLAB_constructor
  end interface

  !!
  !! Printer for ASCII MATLAB output file
  !!
  !! NOTE: Should not check calls logic, which is responsibility of oputputFile class!
  !!
  type, public,extends(asciiOutput) :: asciiMATLAB
    private
    ! State components
    integer(shortInt) :: state  = IN_BLOCK
    type(stackChar)   :: blockNameStack
    integer(shortInt) :: blockLevel = 0

    ! Outputfile
    type(charTape) :: output

    ! Buffers
    type(charTape) :: prefix
  !  character(nameLen)                          :: prefix   = 'R'
    integer(shortInt),dimension(:), allocatable :: shapeBuffer

  contains
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
  !!
  !!
  function asciiMATLAB_constructor() result (new)
    type(asciiMATLAB) :: new
  end function asciiMATLAB_constructor

  !!
  !! Write to a provided unit
  !!
  subroutine writeToFile(self,unit)
    class(asciiMATLAB), intent(inout) :: self
    integer(shortInt), intent(in)     :: unit
    character(:),allocatable          :: form

    form = '(A'//numToChar(self % output % length())//')'

    write(unit,form) self % output % expose()

  end subroutine writeToFile

  !!
  !! Change state to writing new block with "name"
  !!
  subroutine startBlock(self,name)
    class(asciiMATLAB), intent(inout) :: self
    character(nameLen), intent(in)    :: name
    character(100), parameter :: Here ='startBlock (asciiMATLAB_class.f90)'

    ! Check if state support acction
    if (self % state == IN_ENTRY .or. self % state == IN_ARRAY) then
      call fatalError(Here,'Cannot start writing new block inside entry or array')
    end if

    ! Update state - change current prefix
    call self % blockNameStack % push(name)
    call self % prefix % append(trim(name)//'_')
    self % blockLevel = self % blockLevel + 1

  end subroutine startBlock

  !!
  !! End top level block and return to previous block
  !! Return error if this is called in root block
  !!
  subroutine endBlock(self)
    class(asciiMATLAB), intent(inout) :: self
    character(nameLen)                :: temp
    integer(shortInt)                 :: N
    character(100), parameter :: Here ='endBlock (asciiMATLAB_class.f90)'


    ! Check if state support acction
    if (self % state == IN_ENTRY .or. self % state == IN_ARRAY) then
      call fatalError(Here,'Cannot end writing new block inside entry or array')
    end if

    if ( self % blockLevel == 0) then
      call fatalError(Here,'Cannot exit from root block')
    end if

    ! Update state - change current prefix
    call self % blockNameStack % pop (temp)

    ! Calculate how much prefix needs to be cut
    N = len_trim(temp)+1
    call self % prefix % cut(N)
    self % blockLevel = self % blockLevel - 1

  end subroutine endBlock

  !!
  !! Change state to writing a new entry
  !! Can recive single value or array next
  !!
  subroutine startEntry(self,name)
    class(asciiMATLAB), intent(inout) :: self
    character(*), intent(in)          :: name
    character(100), parameter :: Here ='startEntry (asciiMATLAB_class.f90)'

    ! Check if state support acction
    if (self % state == IN_ENTRY .or. self % state == IN_ARRAY) then
      call fatalError(Here,'Cannot star writing new entry inside entry or array')
    end if

    ! Update state
    self % state = IN_ENTRY

    ! Write variable name with prefix
    call self % output % append( trim(self % prefix % expose())//trim(name))
    call self % output % append( ' = ')


  end subroutine startEntry
   !!
  !! End writing a new entry
  !!
  subroutine endEntry(self)
    class(asciiMATLAB), intent(inout) :: self
    character(100), parameter :: Here ='endEntry (asciiMATLAB_class.f90)'

    ! Check if state support acction
    if (self % state == IN_BLOCK .or. self % state == IN_ARRAY) then
      call fatalError(Here,'Cannot finish entry in block or array')
    end if

    ! Update state
    self % state = IN_BLOCK

    ! Write semicolon and newline
    call self % output % append( ';'//NEWLINE)

  end subroutine endEntry

  !!
  !! Start writing array with shape & column-major order(leftmost index varies fastest)
  !! Name should alrady be provided by "startEntry"
  !!
  subroutine startArray(self,shape)
    class(asciiMATLAB), intent(inout)         :: self
    integer(shortInt),dimension(:),intent(in) :: shape
    character(100), parameter :: Here ='startArray (asciiMATLAB_class.f90)'

    ! Check if state support acction
    if (self % state /= IN_ENTRY) then
      call fatalError(Here,'Cannot finish entry in block or array')
    end if

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
  !! End writing array
  !!
  subroutine endArray(self)
    class(asciiMATLAB), intent(inout) :: self
    integer(shortInt)                 :: i
    character(100), parameter :: Here ='endArray (asciiMATLAB_class.f90)'

    ! Check if state support acction
    if (self % state /= IN_ARRAY) then
      call fatalError(Here,'Cannot finish array inside an entry')
    end if

    ! Update state
    self % state = IN_ENTRY

    ! Write end of array
    call self % output % cut(1)
    call self % output % append(']')

    ! Finish reshape function for higher rank arrays
    if( size(self % shapeBuffer) > 1) then
      do i=1,size(self % shapeBuffer)
        call self % output % append(','// numToChar(self % shapeBuffer(i)))
      end do
      call self % output % append(')')
    end if

  end subroutine endArray

  !!
  !! Print val assuming it contains valid printed number with no leading or trailing blanks
  !!
  subroutine printNum(self,val)
    class(asciiMATLAB), intent(inout) :: self
    character(*),intent(in)           :: val
    character(100), parameter :: Here ='printNum (asciiMATLAB_class.f90)'

    if(self % state == IN_ARRAY) then
      call self % output % append(val//',')

    else if( self % state == IN_ENTRY) then
      call self % output % append(val)

    else
      call fatalError(Here,'Cannot print number directly into block')
    end if

  end subroutine printNum

  !!
  !! Print val assuming it contains valid printed string with no leading or trailing blanks
  !!
  subroutine printChar(self,val)
    class(asciiMATLAB), intent(inout) :: self
    character(*),intent(in)           :: val
    character(100), parameter :: Here ='printChar (asciiMATLAB_class.f90)'

    if(self % state == IN_ARRAY) then
      call self % output % append(BRAKET_L//APOS//val//APOS//BRAKET_R //",")

    else if( self % state == IN_ENTRY) then
      call self % output % append(BRAKET_L//APOS//val//APOS//BRAKET_R )

    else
      call fatalError(Here,'Cannot print number directly into block')
    end if

  end subroutine printChar
    
end module asciiMATLAB_class
