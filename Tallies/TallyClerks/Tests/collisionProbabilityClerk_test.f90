module collisionProbabilityClerk_test

  use numPrecision
  use tallyResult_class,               only : tallyResult
  use collisionProbabilityClerk_class, only : collisionProbabilityClerk, CPMResult
  use particle_class,                  only : particle, particleState, P_NEUTRON
  use particleDungeon_class,           only : particleDungeon
  use dictionary_class,                only : dictionary
  use scoreMemory_class,               only : scoreMemory
  use testNeutronDatabase_class,       only : testNeutronDatabase
  use outputFile_class,                only : outputFile
  use funit

  implicit none

@testCase
  type, extends(TestCase) :: test_collisionProbabilityClerk
    private
    type(collisionProbabilityClerk) :: clerk
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_collisionProbabilityClerk

contains

  !!
  !! Sets up test_collisionProbabilityClerk object we can use in a number of tests
  !!
  !! Simple 2x2 collision probability score (with outside region) divided with test map
  !!
  subroutine setUp(this)
    class(test_collisionProbabilityClerk), intent(inout) :: this
    type(dictionary)                                     :: dict
    type(dictionary)                                     :: mapDict
    character(nameLen)                                   :: name

    call mapDict % init(2)
    call mapDict % store('type','testMap')
    call mapDict % store('maxIdx',2)


    ! Build intput dictionary
    call dict % init(2)
    call dict % store('type','collisionProbabilityClerk')
    call dict % store('map', mapDict)

    name = 'testClerk'
    call this % clerk % init(dict,name)


    call mapDict % kill()
    call dict % kill()
  end subroutine setUp

  !!
  !! Kills test_collisionProbabilityClerk object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_collisionProbabilityClerk), intent(inout) :: this

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
    class(test_collisionProbabilityClerk), intent(inout) :: this
    type(scoreMemory)                                    :: mem
    type(particle)                                       :: p
    type(particleDungeon)                                :: pop
    type(testNeutronDatabase)                            :: xsData
    real(defReal)                                        :: val
    class(tallyResult), allocatable                      :: res
    real(defReal), parameter :: TOL = 1.0E-7


    ! Create score memory
    call mem % init(int(this % clerk % getSize(), longInt) , 1, batchSize = 1)
    call this % clerk % setMemAddress(1_longInt)

    ! Create test transport Nuclear Data
    call xsData % build(1.1_defReal, fissionXS = 1.1_defReal, nuFissionXS = 2.0_defReal)

    ! Create one particle that can be made to collide
    ! repeatedly in several materials

    ! Score some events
    p % type = P_NEUTRON

    ! Particle starts in material 2 and collides in material 2
    call p % setMatIdx(2)
    p % w = 0.7
    p % preCollision % matIdx = 2
    call this % clerk % reportInColl(p, xsData, mem)

    ! Particle starts in material 1 and collides in material 2
    call p % setMatIdx(2)
    p % w = 1.1
    p % preCollision % matIdx = 1
    call this % clerk % reportInColl(p, xsData, mem)

    ! Particle starts in material 1 and collides in material 1
    call p % setMatIdx(1)
    p % w = 1.0
    p % preCollision % matIdx = 1
    call this % clerk % reportInColl(p, xsData, mem)

    ! Particle starts in material 2 and collides in material 1
    call p % setMatIdx(1)
    p % w = 1.4
    p % preCollision % matIdx = 2
    call this % clerk % reportInColl(p, xsData, mem)

    ! Particle starts in material 2 and collides in another, unknown material
    call p % setMatIdx(7)
    p % w = 1.0
    p % preCollision % matIdx = 2
    call this % clerk % reportInColl(p, xsData, mem)

    ! Particle starts in an unknown material and collides in material 1
    call p % setMatIdx(1)
    p % w = 0.9
    p % preCollision % matIdx = 88
    call this % clerk % reportInColl(p, xsData, mem)
    call this % clerk % reportCycleEnd(pop, mem)

    ! Close cycle
    call mem % closeCycle(ONE)

    ! Verify results

    ! Collision probability matrix

    ! outside -> outside Transition
    call mem % getResult(val, 1_longInt)
    @assertEqual(ZERO ,val, TOL)

    ! outside -> 1 Transition
    call mem % getResult(val, 2_longInt)
    @assertEqual(ONE ,val, TOL)

    ! outside -> 2 Transition
    call mem % getResult(val, 3_longInt)
    @assertEqual(ZERO ,val, TOL)

    ! 1 -> outside Transition
    call mem % getResult(val, 4_longInt)
    @assertEqual(ZERO ,val, TOL)

    ! 1 -> 1 Transition
    call mem % getResult(val, 5_longInt)
    @assertEqual(0.47619047619_defReal ,val, TOL)

    ! 1 -> 2 Transition
    call mem % getResult(val, 6_longInt)
    @assertEqual(0.52380952381_defReal, val, TOL)

    ! 2 -> outside Transition
    call mem % getResult(val, 7_longInt)
    @assertEqual(0.32258064516_defReal, val, TOL)

    ! 2 -> 1 Transition
    call mem % getResult(val, 8_longInt)
    @assertEqual(0.45161290322_defReal ,val, TOL)

    ! 2 -> 2 Transition
    call mem % getResult(val, 9_longInt)
    @assertEqual(0.22580645161_defReal, val, TOL)

    ! Verify run-time result
    call this % clerk % getResult(res, mem)

    select type(res)
      class is (CPMresult)
        @assertEqual(3, res % N)

        ! outside -> outside Transition
        @assertEqual(ZERO, res  % CPM(1,1,1), TOL)

        ! outside -> 1 Transition
        @assertEqual(ONE, res  % CPM(2,1,1), TOL)

        ! outside -> 2 Transition
        @assertEqual(ZERO, res  % CPM(3,1,1), TOL)

        ! 1 -> outside Transition
        @assertEqual(ZERO, res  % CPM(1,2,1), TOL)

        ! 1 -> 1 Transition
        @assertEqual(0.47619047619_defReal, res  % CPM(2,2,1), TOL)

        ! 1 -> 2 Transition
        @assertEqual(0.52380952381_defReal, res  % CPM(3,2,1), TOL)

        ! 2 -> outside Transition
        @assertEqual(0.32258064516, res  % CPM(1,3,1), TOL)

        ! 2 -> 1 Transition
        @assertEqual(0.45161290322, res  % CPM(2,3,1), TOL)

        ! 2 -> 2 Transition
        @assertEqual(0.22580645161, res  % CPM(3,3,1), TOL)

        ! Clean all entries
        res % CPM = ZERO

      class default
        @assertEqual(1,2)
    end select

    ! Get result again -> verify correcness of reallocation logic by code coverage
    call this % clerk % getResult(res, mem)
    select type(res)
      class is (CPMresult)
        ! 1 -> 1 Transition
        @assertEqual(0.47619047619_defReal, res  % CPM(2,2,1), TOL)

        ! Change size of matrix
        res % N = 2
        deallocate(res % CPM)
        allocate(res % CPM(2,2,1))

    end select

    ! Get result yet again to ensure the size was not incorrectly modified
    call this % clerk % getResult(res, mem)

    select type(res)
      class is (CPMresult)
        @assertEqual(3, res % N)
        @assertEqual([3,3,2], shape(res % CPM))
    end select

    ! Clean
    call xsData % kill()
    call pop % kill()

  end subroutine testSimpleUseCase



  !!
  !! Test correctness of the printing calls
  !!
@Test
  subroutine testPrintingCorrectness(this)
    class(test_collisionProbabilityClerk), intent(inout) :: this
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

end module collisionProbabilityClerk_test
