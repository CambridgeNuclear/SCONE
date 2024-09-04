module shannonEntropyClerk_test

  use numPrecision
  use shannonEntropyClerk_class,      only : shannonEntropyClerk
  use particle_class,                 only : particle, particleState
  use particleDungeon_class,          only : particleDungeon
  use dictionary_class,               only : dictionary
  use scoreMemory_class,              only : scoreMemory
  use outputFile_class,               only : outputFile
  use funit

  implicit none

@testCase
  type, extends(TestCase) :: test_shannonEntropyClerk
    private
    type(shannonEntropyClerk) :: clerk
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_shannonEntropyClerk

contains

  !!
  !! Sets up test_simpleFMClerk object we can use in a number of tests
  !!
  !! Simple 3x3 fission matrix divided with test map
  !!
  subroutine setUp(this)
    class(test_shannonEntropyClerk), intent(inout) :: this
    type(dictionary)                               :: dict
    type(dictionary)                               :: mapDict
    character(nameLen)                             :: name

    call mapDict % init(2)
    call mapDict % store('type','testMap')
    call mapDict % store('maxIdx',2)


    ! Build intput dictionary
    call dict % init(2)
    call dict % store('type','shannonEntropyClerk')
    call dict % store('map', mapDict)
    call dict % store('cycles', 2)

    name = 'testClerk'
    call this % clerk % init(dict,name)

    call mapDict % kill()
    call dict % kill()
  end subroutine setUp

  !!
  !! Kills test_shannonEntropyClerk object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_shannonEntropyClerk), intent(inout) :: this

    call this % clerk % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>


  !!
  !! Test correctness in a simple use case
  !!
@Test
  subroutine testSimpleUseCase(this)
    class(test_shannonEntropyClerk), intent(inout) :: this
    type(scoreMemory)                              :: mem
    type(particleState)                            :: phase
    type(particleDungeon)                          :: pop
    real(defReal), parameter :: TOL = 1.0E-7

    ! Create score memory
    call mem % init(int(this % clerk % getSize(), longInt) , 1, batchSize = 1)
    call this % clerk % setMemAddress(1_longInt)

    ! Crate dungeon of original events
    ! One particle born in matIdx 1 and other in 2
    call pop % init(3)

    phase % wgt = ONE
    phase % matIdx = 2
    call pop % detain(phase)

    phase % wgt = ONE
    phase % matIdx = 1
    call pop % detain(phase)

    call this % clerk % reportCycleEnd(pop, mem)

    ! Close cycle
    call mem % closeCycle(ONE)

    ! Verify results for uniform distribution
    @assertEqual(ONE, this % clerk % value(1), TOL)

    ! Move particles to the same bin
    call pop % kill()
    call pop % init(3)
    call pop % detain(phase)
    call pop % detain(phase)

    call this % clerk % reportCycleEnd(pop,  mem)

    ! Close cycle
    call mem % closeCycle(ONE)

    ! Verify results for all particles in one bine
    @assertEqual(ZERO, this % clerk % value(2), TOL)

    ! Clean
    call pop % kill()
  end subroutine testSimpleUseCase



  !!
  !! Test correctness of the printing calls
  !!
@Test
  subroutine testPrintingCorrectness(this)
    class(test_shannonEntropyClerk), intent(inout) :: this
    type(outputFile)                               :: outF
    type(scoreMemory)                              :: mem

    ! Create score memory
    call mem % init(int(this % clerk % getSize(), longInt) , 1)
    call this % clerk % setMemAddress(1_longInt)

    ! Verify that output calls are correct
    call outF % init('dummyPrinter', fatalErrors = .false.)
    call this % clerk % print (outF, mem)

    @assertTrue(outF % isValid())

  end subroutine testPrintingCorrectness

end module shannonEntropyClerk_test
