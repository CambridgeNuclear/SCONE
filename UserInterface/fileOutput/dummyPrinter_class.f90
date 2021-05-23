module dummyPrinter_class

  use numPrecision
  use genericProcedures, only : fatalError, numToChar
  use stack_class,       only : stackChar
  use asciiOutput_inter, only : asciiOutput
  use charTape_class,    only : charTape

  implicit none
  private

  !!
  !! Dummy printer
  !!
  !! Accepts printer calls but does not produce output.
  !! Used for testing.
  !!
  !! NOTE:
  !!  Does not check calls sequence, which is responsibility of oputputFile class!
  !!
  type, public,extends(asciiOutput) :: dummyPrinter
    private

  contains
    procedure :: init
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
  !! Initialise the printer
  !!
  !! See asciiOutput_inter for details
  !!
  subroutine init(self)
    class(dummyPrinter), intent(inout) :: self

    ! Nothing to do

  end subroutine init

  !!
  !! Print the output to the given unit
  !!
  !! See asciiOutput_inter for details
  !!
  subroutine writeToFile(self, unit)
    class(dummyPrinter), intent(inout) :: self
    integer(shortInt), intent(in)      :: unit


  end subroutine writeToFile

  !!
  !! Change state to writing new block with "name"
  !!
  !! See asciiOutput_inter for details
  !!
  subroutine startBlock(self, name)
    class(dummyPrinter), intent(inout) :: self
    character(nameLen), intent(in)    :: name

  end subroutine startBlock

  !!
  !! End current block
  !!
  !! See asciiOutput_inter for details
  !!
  subroutine endBlock(self)
    class(dummyPrinter), intent(inout) :: self


  end subroutine endBlock

  !!
  !! Change state to writing a new entry
  !!
  !! See asciiOutput_inter for details
  !!
  subroutine startEntry(self, name)
    class(dummyPrinter), intent(inout) :: self
    character(*), intent(in)           :: name


  end subroutine startEntry

  !!
  !! End writing a new entry
  !!
  !! See asciiOutput_inter for details
  !!
  subroutine endEntry(self)
    class(dummyPrinter), intent(inout) :: self

  end subroutine endEntry

  !!
  !! Start writing array with the given shape
  !!
  !! See asciiOutput_inter for details
  !!
  subroutine startArray(self, shape)
    class(dummyPrinter), intent(inout) :: self
    integer(shortInt),dimension(:),intent(in)     :: shape

  end subroutine startArray

  !!
  !! End writing the array
  !!
  !! See asciiOutput_inter for details
  !!
  subroutine endArray(self)
    class(dummyPrinter), intent(inout) :: self


  end subroutine endArray

  !!
  !! Print numerical value
  !!
  !! See asciiOutput_inter for details
  !!
  subroutine printNum(self, val)
    class(dummyPrinter), intent(inout) :: self
    character(*),intent(in)           :: val


  end subroutine printNum

  !!
  !! Print character value
  !!
  !! See asciiOutput_inter for details
  !!
  subroutine printChar(self, val)
    class(dummyPrinter), intent(inout) :: self
    character(*),intent(in)           :: val

  end subroutine printChar

end module dummyPrinter_class
