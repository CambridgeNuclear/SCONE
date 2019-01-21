module dummyPrinter_class

  use numPrecision
  use genericProcedures, only : fatalError, numToChar
  use stack_class,       only : stackChar
  use asciiOutput_inter, only : asciiOutput
  use charTape_class,    only : charTape

  implicit none
  private

  !!
  !! Constructor
  !!
  interface dummyPrinter
    module procedure dummyPrinter_constructor
  end interface

  !!
  !! Printer for ASCII MATLAB output file
  !!
  !! NOTE: Should not check calls logic, which is responsibility of oputputFile class!
  !!
  type, public,extends(asciiOutput) :: dummyPrinter
    private

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

  end type dummyPrinter

contains

  !!
  !!
  !!
  function dummyPrinter_constructor() result (new)
    type(dummyPrinter) :: new
  end function dummyPrinter_constructor

  !!
  !! Write to a provided unit
  !!
  subroutine writeToFile(self,unit)
    class(dummyPrinter), intent(inout) :: self
    integer(shortInt), intent(in)      :: unit


  end subroutine writeToFile

  !!
  !! Change state to writing new block with "name"
  !!
  subroutine startBlock(self,name)
    class(dummyPrinter), intent(inout) :: self
    character(nameLen), intent(in)    :: name

  end subroutine startBlock

  !!
  !! End top level block and return to previous block
  !! Return error if this is called in root block
  !!
  subroutine endBlock(self)
    class(dummyPrinter), intent(inout) :: self


  end subroutine endBlock

  !!
  !! Change state to writing a new entry
  !! Can recive single value or array next
  !!
  subroutine startEntry(self,name)
    class(dummyPrinter), intent(inout) :: self
    character(*), intent(in)           :: name


  end subroutine startEntry

  !!
  !! End writing a new entry
  !!
  subroutine endEntry(self)
    class(dummyPrinter), intent(inout) :: self

  end subroutine endEntry

  !!
  !! Start writing array with shape & column-major order(leftmost index varies fastest)
  !! Name should alrady be provided by "startEntry"
  !!
  subroutine startArray(self,shape)
    class(dummyPrinter), intent(inout) :: self
    integer(shortInt),dimension(:),intent(in)     :: shape

  end subroutine startArray

  !!
  !! End writing array
  !!
  subroutine endArray(self)
    class(dummyPrinter), intent(inout) :: self


  end subroutine endArray

  !!
  !! Print val assuming it contains valid printed number with no leading or trailing blanks
  !!
  subroutine printNum(self,val)
    class(dummyPrinter), intent(inout) :: self
    character(*),intent(in)           :: val


  end subroutine printNum

  !!
  !! Print val assuming it contains valid printed string with no leading or trailing blanks
  !!
  subroutine printChar(self,val)
    class(dummyPrinter), intent(inout) :: self
    character(*),intent(in)           :: val

  end subroutine printChar
    
end module dummyPrinter_class
