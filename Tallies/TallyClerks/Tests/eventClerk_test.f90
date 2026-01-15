module eventClerk_test

  use numPrecision
  use eventClerk_class,          only : eventClerk
  use particle_class,            only : particle, particleState, P_NEUTRON
  use dictionary_class,          only : dictionary
  use scoreMemory_class,         only : scoreMemory
  use testNeutronDatabase_class, only : testNeutronDatabase
  use outputFile_class,          only : outputFile
  use funit

  implicit none

@testCase
  type, extends(TestCase) :: test_eventClerk
    private
    type(eventClerk) :: clerk
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_eventClerk

contains

  !!
  !! Sets up test_eventClerk object we can use in a number of tests
  !!
  !! Simple 2x2 collision probability score (with outside region) divided with test map
  !!
  subroutine setUp(this)
    class(test_eventClerk), intent(inout) :: this
    type(dictionary)                      :: dict
    type(dictionary)                      :: mapDict
    character(nameLen)                    :: name

    call mapDict % init(2)
    call mapDict % store('type','testMap')
    call mapDict % store('maxIdx',1)

    ! Build intput dictionary
    call dict % init(3)
    call dict % store('type','eventClerk')
    call dict % store('file','file.txt')
    call dict % store('map', mapDict)

    name = 'testClerk'
    call this % clerk % init(dict,name)

    call mapDict % kill()
    call dict % kill()
  end subroutine setUp

  !!
  !! Kills test_eventClerk object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_eventClerk), intent(inout) :: this

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
    class(test_eventClerk), intent(inout) :: this
    type(scoreMemory)                     :: mem
    type(particle)                        :: p
    type(testNeutronDatabase)             :: xsData
    real(defReal)                         :: valReal
    real(defReal), dimension(3)           :: arrReal
    integer(shortInt)                     :: valInt
    integer(shortInt)                     :: MT
    real(defReal), parameter :: TOL = 1.0E-5

    ! Create score memory
    ! Only a dummy since it is not directly used by eventClerk
    call mem % init(1_longInt , 1, batchSize = 1)
    call this % clerk % setMemAddress(1_longInt)

    ! Create test transport Nuclear Data
    call xsData % build(1.1_defReal, fissionXS = 1.1_defReal, nuFissionXS = 2.0_defReal)

    ! Create one particle that can be made to collide repeatedly

    ! Score some events
    p % type = P_NEUTRON

    ! Particle collides in material 1 and loses energy
    p % w = 0.7
    p % broodID = 1_shortInt
    p % time = -ONE
    call p % point([ZERO, ZERO, ONE])
    call p % teleport([-3.0_defReal, 10.0_defReal, 7.2_defReal])
    p % preCollision % E = 7.0_defReal
    p % E = 5.0_defReal
    MT = 53
    p % isDead = .false.
    call p % setMatIdx(1)
    call this % clerk % reportOutColl(p, MT, ONE, xsData, mem)

    ! Particle collides in material 2: this event shouldn't register!
    p % w = 1.1
    p % broodID = 5_shortInt
    p % time = 10.0_defReal
    p % preCollision % E = 70.0_defReal
    p % E = 1.0_defReal
    MT = 7
    p % isDead = .false.
    call p % setMatIdx(2)
    call this % clerk % reportOutColl(p, MT, ONE, xsData, mem)

    ! Back in material 1 with a different set of parameters
    p % w = 4.0
    p % broodID = 11_shortInt
    p % time = 2.2_defReal
    call p % point([ZERO, ONE, ZERO])
    call p % teleport([5.2_defReal, -1.0_defReal, -9.0_defReal])
    p % preCollision % E = 0.1_defReal
    p % E = 1.8_defReal
    MT = -8
    p % isDead = .false.
    call p % setMatIdx(1)
    call this % clerk % reportOutColl(p, MT, ONE, xsData, mem)

    ! Verify results
    ! First collision
    valReal = real(this % clerk % scores(1,1) % Eincident, defReal)
    @assertEqual(7.0_defReal, valReal, TOL)
    valReal = real(this % clerk % scores(1,1) % Edeposit, defReal)
    @assertEqual(2.0_defReal, valReal, TOL)
    valReal = real(this % clerk % scores(1,1) % time, defReal)
    @assertEqual(-1.0_defReal, valReal, TOL)
    valReal = real(this % clerk % scores(1,1) % w, defReal)
    @assertEqual(0.7_defReal, valReal, TOL)
    arrReal = real(this % clerk % scores(1,1) % r, defReal)
    @assertEqual([-3.0_defReal, 10.0_defReal, 7.2_defReal], arrReal, TOL)
    arrReal = real(this % clerk % scores(1,1) % dir, defReal)
    @assertEqual([ZERO, ZERO, ONE], arrReal, TOL)
    valInt = this % clerk % scores(1,1) % brood
    @assertEqual(1, valInt)
    valInt = this % clerk % scores(1,1) % MT
    @assertEqual(53, valInt)

    ! Second collision
    valReal = real(this % clerk % scores(2,1) % time, defReal)
    @assertEqual(2.2_defReal, valReal, TOL)
    valReal = real(this % clerk % scores(2,1) % Eincident, defReal)
    @assertEqual(0.1_defReal, valReal, TOL)
    valReal = real(this % clerk % scores(2,1) % Edeposit, defReal)
    @assertEqual(-1.7_defReal, valReal, TOL)
    valReal = real(this % clerk % scores(2,1) % w, defReal)
    @assertEqual(4.0_defReal, valReal, TOL)
    arrReal = real(this % clerk % scores(2,1) % r, defReal)
    @assertEqual([5.2_defReal, -1.0_defReal, -9.0_defReal], arrReal, TOL)
    arrReal = real(this % clerk % scores(2,1) % dir, defReal)
    @assertEqual([ZERO, ONE, ZERO], arrReal, TOL)
    valInt = this % clerk % scores(2,1) % brood
    @assertEqual(11, valInt)
    valInt = this % clerk % scores(2,1) % MT
    @assertEqual(-8, valInt)

    ! Clean
    call xsData % kill()

  end subroutine testSimpleUseCase

  !!
  !! Test correctness of the printing calls
  !!
@Test
  subroutine testPrintingCorrectness(this)
    class(test_eventClerk), intent(inout) :: this
    type(outputFile)                     :: outF
    type(scoreMemory)                    :: mem

    ! Create score memory
    call mem % init(1_longInt , 1)
    call this % clerk % setMemAddress(1_longInt)

    ! Verify that output calls are correct
    call outF % init('dummyPrinter', fatalErrors = .false.)
    call this % clerk % print(outF, mem)

    @assertTrue(outF % isValid())

  end subroutine testPrintingCorrectness

end module eventClerk_test
