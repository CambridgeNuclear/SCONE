module simpleFMClerk_test

  use numPrecision
  use tallyResult_class,              only : tallyResult
  use simpleFMClerk_class,            only : simpleFMClerk
  use particle_class,                 only : particle, particleState
  use particleDungeon_class,          only : particleDungeon
  use dictionary_class,               only : dictionary
  use scoreMemory_class,              only : scoreMemory
  use testTransportNuclearData_class, only : testTransportNuclearData
  use outputFile_class,               only : outputFile
  use pFUnit_mod

  implicit none

@testCase
  type, extends(TestCase) :: test_simpleFMClerk
    private
    type(simpleFMClerk) :: clerk
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_simpleFMClerk

contains

  !!
  !! Sets up test_simpleFMClerk object we can use in a number of tests
  !!
  !! Simple 3x3 fission matrix divided with test map
  !!
  subroutine setUp(this)
    class(test_simpleFMClerk), intent(inout) :: this
    type(dictionary)                         :: dict
    type(dictionary)                         :: mapDict
    character(nameLen)                       :: name

    call mapDict % init(2)
    call mapDict % store('type','testMap')
    call mapDict % store('maxIdx',3)


    ! Build intput dictionary
    call dict % init(2)
    call dict % store('type','simpleFMClerk')
    call dict % store('map', mapDict)

    name = 'testClerk'
    call this % clerk % init(dict,name)


    call mapDict % kill()
    call dict % kill()
  end subroutine setUp

  !!
  !! Kills test_simpleFMClerk object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_simpleFMClerk), intent(inout) :: this

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
    class(test_simpleFMClerk), intent(inout) :: this
    type(scoreMemory)                        :: mem
    type(particle)                           :: p
    type(particleState)                      :: phase
    type(particleDungeon)                    :: pop
    type(testTransportNuclearData),pointer   :: xsData
    real(defReal)                            :: val
    real(defReal), parameter :: TOL = 1.0E-7

    ! Create score memory
    call mem % init(int(this % clerk % getSize(), longInt) , 1, batchSize = 1)
    call this % clerk % setMemAddress(1_longInt)

    ! Create test transport Nuclear Data
    allocate(xsData)
    call xsData % build(1.1_defReal, fissionXS = 1.1_defReal, nuFissionXS = 2.0_defReal)
    p % xsData => xsData

    ! Crate dungeon of original events
    ! One particle born in matIdx 1 and other in 2
    call pop % init(3)

    phase % wgt = ONE
    phase % matIdx = 2
    call pop % detain(phase)

    phase % wgt = ONE
    phase % matIdx = 1
    call pop % detain(phase)

    call this % clerk % reportCycleStart(pop, mem)

    ! Score some events

    call p % setMatIdx(2)
    p % w = 0.7
    p % preHistory % matIdx = 2
    call this % clerk % reportInColl(p, mem)

    call p % setMatIdx(1)
    p % w = 1.1
    p % preHistory % matIdx = 2
    call this % clerk % reportInColl(p, mem)

    call p % setMatIdx(1)
    p % w = 1.0
    p % preHistory % matIdx = 1
    call this % clerk % reportInColl(p, mem)

    call this % clerk % reportCycleEnd(pop, mem)

    ! Close cycle
    call mem % closeCycle(ONE)

    ! Verify results

    ! Fission matrix
    ! 1 -> 1 Transition
    call mem % getResult(val, 1_longInt)
    @assertEqual(1.818181818181_defReal ,val, TOL)

    ! 1 -> 2 Transition
    call mem % getResult(val, 2_longInt)
    @assertEqual(ZERO, val, TOL)

    ! 1 -> 3 Transition
    call mem % getResult(val, 3_longInt)
    @assertEqual(ZERO, val, TOL)

    ! 2 -> 1 Transition
    call mem % getResult(val, 4_longInt)
    @assertEqual(2.0_defReal, val, TOL)

    ! 2 -> 2 Transition
    call mem % getResult(val, 5_longInt)
    @assertEqual(1.27272727272727_defReal, val, TOL)

    ! Clean
    call xsData % kill()
    call pop % kill()
    deallocate(xsData)
  end subroutine testSimpleUseCase



  !!
  !! Test correctness of the printing calls
  !!
@Test
  subroutine testPrintingCorrectness(this)
    class(test_simpleFMClerk), intent(inout) :: this
    type(outputFile)                         :: outF
    type(scoreMemory)                        :: mem

    ! Create score memory
    call mem % init(int(this % clerk % getSize(), longInt) , 1)
    call this % clerk % setMemAddress(1_longInt)

    ! Verify that output calls are correct
    call outF % init('dummyPrinter', fatalErrors = .false.)
    call this % clerk % print (outF, mem)

    @assertTrue(outF % isValid())

  end subroutine testPrintingCorrectness



end module simpleFMClerk_test
