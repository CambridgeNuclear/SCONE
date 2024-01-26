module outputFile_test

  use numPrecision
  use outputFile_class, only : outputFile
  use funit

  implicit none

@testCase
  type, extends(TestCase) :: test_outputFile
    private
    type(outputFile) :: outFile
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_outputFile


contains

  !!
  !! Sets up test_outputFile object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_outputFile), intent(inout) :: this

    call this % outFile % init('dummyPrinter', fatalErrors = .false.)

  end subroutine setUp

  !!
  !! Kills test_outputFile object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_outputFile), intent(inout) :: this

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test that block errors are caught
  !!
@Test
  subroutine testBlockLogic(this)
    class(test_outputFile), intent(inout) :: this
    character(nameLen)                    :: myBlock, myBlock2, nextBlock, myArray

    myBlock = 'myBlock'
    myBlock2 = 'myBlock2'
    nextBlock = 'nextBlock'
    myArray = 'myArray'

    ! Test for too early block closure
    call this % outFile % endBlock()
    @assertFalse( this % outFile % isValid())
    call this % outFile % reset()

    ! Test for to many block closures
    call this % outFile % startBlock(myBlock)
    call this % outFile % startBlock(myBlock2)
    call this % outFile % endBlock()
    call this % outFile % endBlock()
    call this % outFile % endBlock()
    @assertFalse( this % outFile % isValid())
    call this % outFile % reset()

    ! Try to close block from an array
    call this % outFile % startBlock(myBlock)
    call this % outFile % startArray(myArray,[2,3])
    call this % outFile % endBlock()
    @assertFalse( this % outFile % isValid())
    call this % outFile % reset()

    ! Try to start block from an array
    call this % outFile % startBlock(myBlock)
    call this % outFile % startArray(myArray,[2,3])
    call this % outFile % startBlock(nextBlock)
    @assertFalse( this % outFile % isValid())
    call this % outFile % reset()

  end subroutine testBlockLogic

  !!
  !! Test that Array Errors are caught
  !!
@Test
  subroutine testArrayLogic(this)
    class(test_outputFile), intent(inout) :: this
    character(nameLen), dimension(2) :: charArray
    character(nameLen)               :: name
    character(nameLen)               :: myBlock, myArray, myArray2
    integer(shortInt), dimension(:), allocatable :: temp_int

    myBlock = 'myBlock'
    myArray = 'myArray'
    myArray2 = 'myArray2'

    ! Start array in array
    call this % outFile % startBlock(myBlock)
    call this % outFile % startArray(myArray,[2,3])
    call this % outFile % startArray(myArray2,[2])
    @assertFalse( this % outFile % isValid())
    call this % outFile % reset()

    ! Change from result to value
    call this % outFile % startBlock(myBlock)
    call this % outFile % startArray(myArray,[3])
    call this % outFile % addResult(0.0_defReal, 1.0_defReal)
    call this % outFile % addValue(0.0_defReal)
    @assertFalse( this % outFile % isValid())
    call this % outFile % reset()

    ! Change type value in value array Real -> Int
    call this % outFile % startBlock(myBlock)
    call this % outFile % startArray(myArray,[3])
    call this % outFile % addValue(0.0_defReal)
    call this % outFile % addValue(0)
    @assertFalse( this % outFile % isValid())
    call this % outFile % reset()

    ! Change type value in value array Real -> Char
    call this % outFile % startBlock(myBlock)
    call this % outFile % startArray(myArray,[3])
    call this % outFile % addValue(0.0_defReal)
    name = 'jgvjj'
    call this % outFile % addValue(name)
    @assertFalse( this % outFile % isValid())
    call this % outFile % reset()

    ! Change type value in value array Int -> Real
    call this % outFile % startBlock(myBlock)
    call this % outFile % startArray(myArray,[3])
    call this % outFile % addValue(0)
    call this % outFile % addValue(0.0_defReal)
    @assertFalse( this % outFile % isValid())
    call this % outFile % reset()

    ! Array overflow defReal
    call this % outFile % startBlock(myBlock)
    call this % outFile % startArray(myArray,[2])
    call this % outFile % addValue(0.0_defReal)
    call this % outFile % addValue(0.0_defReal)
    call this % outFile % addValue(0.0_defReal)
    @assertFalse( this % outFile % isValid())
    call this % outFile % reset()

    ! Array overflow shortInt
    call this % outFile % startBlock(myBlock)
    call this % outFile % startArray(myArray,[2])
    call this % outFile % addValue(0)
    call this % outFile % addValue(0)
    call this % outFile % addValue(0)
    @assertFalse( this % outFile % isValid())
    call this % outFile % reset()

    ! Array overflow longInt
    call this % outFile % startBlock(myBlock)
    call this % outFile % startArray(myArray,[2])
    call this % outFile % addValue(0_longInt)
    call this % outFile % addValue(0_longInt)
    call this % outFile % addValue(0_longInt)
    @assertFalse( this % outFile % isValid())
    call this % outFile % reset()

    ! Array overflow char
    name ='lkm'
    call this % outFile % startBlock(myBlock)
    call this % outFile % startArray(myArray,[2])
    call this % outFile % addValue(name)
    call this % outFile % addValue(name)
    call this % outFile % addValue(name)
    @assertFalse( this % outFile % isValid())
    call this % outFile % reset()

    ! Array overflow Result
    call this % outFile % startBlock(myBlock)
    call this % outFile % startArray(myArray,[2])
    call this % outFile % addResult(0.0_defReal, 1.0_defReal)
    call this % outFile % addResult(0.0_defReal, 1.0_defReal)
    call this % outFile % addResult(0.0_defReal, 1.0_defReal)
    @assertFalse( this % outFile % isValid())
    call this % outFile % reset()

    ! Close array with too little entries provided
    call this % outFile % startBlock(myBlock)
    call this % outFile % startArray(myArray,[2])
    call this % outFile % addValue(0.0_defReal)
    call this % outFile % endArray()
    @assertFalse( this % outFile % isValid())
    call this % outFile % reset()

    ! Printing shortInt value inside an array
    call this % outFile % startBlock(myBlock)
    call this % outFile % startArray(myArray,[2])
    call this % outFile % addValue(0.0_defReal)
    name ='myVal'
    call this % outFile % printValue(0, name)
    @assertFalse( this % outFile % isValid())
    call this % outFile % reset()

    ! Printing longInt value inside an array
    call this % outFile % startBlock(myBlock)
    call this % outFile % startArray(myArray,[2])
    call this % outFile % addValue(0.0_defReal)
    name ='myVal'
    call this % outFile % printValue(0_longInt, name)
    @assertFalse( this % outFile % isValid())
    call this % outFile % reset()

    ! Printing defReal value inside an array
    call this % outFile % startBlock(myBlock)
    call this % outFile % startArray(myArray,[2])
    call this % outFile % addValue(0.0_defReal)
    name ='myVal'
    call this % outFile % printValue(0.0_defReal, name)
    @assertFalse( this % outFile % isValid())
    call this % outFile % reset()

    ! Printing character value inside an array
    charArray(1) = 'sth'
    charArray(2) ='sth else'
    call this % outFile % startBlock(myBlock)
    call this % outFile % startArray(myArray,[2])
    call this % outFile % addValue(0.0_defReal)
    name ='myVal'
    call this % outFile % printValue(charArray(1), name)
    @assertFalse( this % outFile % isValid())
    call this % outFile % reset()

    ! Printing result inside an array
    call this % outFile % startBlock(myBlock)
    call this % outFile % startArray(myArray,[2])
    call this % outFile % addValue(0.0_defReal)
    name ='myVal'
    call this % outFile % printResult(0.0_defReal,1.0_defReal, name)
    @assertFalse( this % outFile % isValid())
    call this % outFile % reset()

    ! Add defReal value without starting an array
    call this % outFile % startBlock(myBlock)
    call this % outFile % addValue(0.0_defReal)
    @assertFalse( this % outFile % isValid())
    call this % outFile % reset()

    ! Add shortInt value without starting an array
    call this % outFile % startBlock(myBlock)
    call this % outFile % addValue(0)
    @assertFalse( this % outFile % isValid())
    call this % outFile % reset()

    ! Add longInt value without starting an array
    call this % outFile % startBlock(myBlock)
    call this % outFile % addValue(0_longInt)
    @assertFalse( this % outFile % isValid())
    call this % outFile % reset()

    ! Add char value without starting an array
    name ='char value'
    call this % outFile % startBlock(myBlock)
    call this % outFile % addValue(name)
    @assertFalse( this % outFile % isValid())
    call this % outFile % reset()

    ! Add result without starting an array
    call this % outFile % startBlock(myBlock)
    call this % outFile % addResult(0.0_defReal, 0.0_defReal)
    @assertFalse( this % outFile % isValid())
    call this % outFile % reset()

    ! Create a degenerate arrays
    allocate(temp_int(0))
    call this % outFile % startArray(name, temp_int)
    call this % outFile % endArray()
    @assertFalse( this % outFile % isValid())
    call this % outFile % reset()

    call this % outFile % startArray(name, [2, 0, 7])
    call this % outFile % endArray()
    @assertFalse( this % outFile % isValid())
    call this % outFile % reset()

  end subroutine testArrayLogic

  !!
  !! Test Repeated Names
  !!
@Test
  subroutine testRepeatedNames(this)
    class(test_outputFile), intent(inout) :: this
    character(nameLen)                    :: name

    name = 'myKey'

    ! Print ordinary values
    call this % outFile % printResult(1.0_defReal, 0.5_defReal, name)
    call this % outFile % printResult(1.0_defReal, 1.5_defReal, name)
    @assertFalse( this % outFile % isValid())
    call this % outFile % reset()

    call this % outFile % printValue(1.0_defReal, name)
    call this % outFile % printValue(2.0_defReal, name)
    @assertFalse( this % outFile % isValid())
    call this % outFile % reset()

    call this % outFile % printValue(1, name)
    call this % outFile % printValue(2, name)
    @assertFalse( this % outFile % isValid())
    call this % outFile % reset()

    call this % outFile % printValue(1_longInt, name)
    call this % outFile % printValue(2_longInt, name)
    @assertFalse( this % outFile % isValid())
    call this % outFile % reset()

    call this % outFile % printValue("char", name)
    call this % outFile % printValue("charizard", name)
    @assertFalse( this % outFile % isValid())
    call this % outFile % reset()

    ! Test with blocks & arrays
    call this % outFile % startBlock(name)
    call this % outFile % endBlock()
    call this % outFile % startBlock(name)
    call this % outFile % endBlock()
    @assertFalse( this % outFile % isValid())
    call this % outFile % reset()

    call this % outFile % startBlock(name)
    call this % outFile % endBlock()
    call this % outFile % startArray(name, [1])
    call this % outFile % addValue(1)
    call this % outFile % endArray()
    @assertFalse( this % outFile % isValid())
    call this % outFile % reset()

  end subroutine testRepeatedNames

  !!
  !! Test deep nested blocks
  !!
@Test
  subroutine testNestedBlocks(this)
    class(test_outputFile), intent(inout) :: this
    character(nameLen)                    :: name
    integer(shortInt)                     :: i

    name = 'myKey'

    do i = 1, 30
      call this % outFile % startBlock(name)
    end do
    do i = 1, 30
      call this % outFile % endBlock()
    end do

    @assertTrue( this % outFile % isValid())
    call this % outFile % reset()

  end subroutine testNestedBlocks


end module outputFile_test
