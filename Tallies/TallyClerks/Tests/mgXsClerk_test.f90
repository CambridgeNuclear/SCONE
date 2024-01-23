module mgXsClerk_test

  use numPrecision
  use endfConstants
  use genericProcedures,         only : numToChar
  use mgXsClerk_class,           only : mgXsClerk
  use particle_class,            only : particle, particleState
  use dictionary_class,          only : dictionary
  use scoreMemory_class,         only : scoreMemory
  use testNeutronDatabase_class, only : testNeutronDatabase
  use outputFile_class,          only : outputFile
  use pFUnit_mod

  implicit none

  @testCase
    type, extends(TestCase) :: test_mgXsClerk
      private
      type(mgXsClerk)  :: clerk_test1
      type(mgXsClerk)  :: clerk_test2
      type(testNeutronDatabase) :: nucData
    contains
      procedure :: setUp
      procedure :: tearDown
  end type test_mgXsClerk

contains

  !!
  !! Sets up test_mgXsClerk object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_mgXsClerk), intent(inout) :: this
    type(dictionary)                     :: tempDict, energyDict, spaceDict
    character(nameLen)                   :: temp

    ! Build energy map
    call energyDict % init(3)
    call energyDict % store('type', 'energyMap')
    call energyDict % store('grid', 'unstruct')
    call energyDict % store('bins', [1.0E-03_defReal, ONE, 10.0_defReal])

    ! Build material map
    call spaceDict % init(2)
    call spaceDict % store('type', 'testMap')
    call spaceDict % store('maxIdx', 2)

    ! Define first clerk, with high order scattering and spatial map
    call tempDict % init(2)
    call tempDict % store('energyMap', energyDict)
    call tempDict % store('spaceMap', spaceDict)

    temp = 'MGxs1'
    call this % clerk_test1 % init(tempDict, temp)
    call spaceDict % kill()
    call energyDict % kill()
    call tempDict % kill()

    ! Build energy map
    call energyDict % init(3)
    call energyDict % store('type', 'energyMap')
    call energyDict % store('grid', 'unstruct')
    call energyDict % store('bins', [1.0E-11_defReal, 0.6_defReal, 1.2_defReal, 20.0_defReal])

    ! Define second clerk, without spatial map and high order scattering
    call tempDict % init(2)
    call tempDict % store('energyMap', energyDict)
    call tempDict % store('PN', 0)

    temp = 'MGxs2'
    call this % clerk_test2 % init(tempDict, temp)
    call energyDict % kill()
    call tempDict % kill()

    ! Build test neutronDatabase
    call this % nucData % build(ONE, captureXS = 2.0_defReal, &
                                fissionXS = 1.5_defReal, nuFissionXS = 3.0_defReal)

  end subroutine setUp

  !!
  !! Kills test case object
  !!
  subroutine tearDown(this)
    class(test_mgXsClerk), intent(inout) :: this

    call this % clerk_test1 % kill()
    call this % clerk_test2 % kill()
    call this % nucData % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Scoring test for clerk 1
  !!
@Test
  subroutine testScoring_clerk1(this)
    class(test_mgXsClerk), intent(inout) :: this
    character(:),allocatable             :: case
    type(scoreMemory)                    :: mem
    type(particle)                       :: p
    type(particleState)                  :: pFiss
    type(outputFile)                     :: out
    real(defReal), dimension(:,:), allocatable :: fiss, capt, transFL, transOS, &
                                                  nu, chi, P0, P1, P2, P3, P4,  &
                                                  P5, P6, P7, prod
    real(defReal), parameter :: TOL = 1.0E-9

    ! Configure memory
    call mem % init(1000_longInt, 1)
    call this % clerk_test1 % setMemAddress(1_longInt)

    ! Report fission particles
    p % isMG = .false.
    p % w    = 0.5_defReal
    p % E    = 0.3_defReal
    call p % setMatIdx(2)

    ! Scoring
    pFiss = p
    call this % clerk_test1 % reportSpawn(p, pFiss, this % nucData, mem)

    p % E = 3.0_defReal

    ! Scoring
    pFiss = p
    call this % clerk_test1 % reportSpawn(p, pFiss, this % nucData, mem)

    ! Scoring
    call this % clerk_test1 % reportInColl(p, this % nucData, mem, .false.)

    p % preCollision % wgt = 0.2_defReal
    p % preCollision % E   = 3.0_defReal
    p % preCollision % matIdx = 2
    p % preCollision % isMG   = .false.
    p % E = 0.1_defReal

    call this % clerk_test1 % reportOutColl(p, N_2N, 0.75_defReal, this % nucData, mem)
    call mem % closeCycle(ONE)

    ! Process and get results
    call this % clerk_test1 % processRes(mem, capt, fiss, transFL, transOS, nu, chi, P0, P1, prod)
    call this % clerk_test1 % processPN(mem, P2, P3, P4, P5, P6, P7)

    ! Verify results of scoring
    @assertEqual([ZERO, ZERO, TWO, ZERO], capt(1,:), TOL, 'Capture XS' )
    @assertEqual([ZERO, ZERO, 1.5_defReal, ZERO], fiss(1,:), TOL, 'Fission XS' )
    @assertEqual([ZERO, ZERO, TWO, ZERO], nu(1,:), TOL, 'NuFission XS' )
    @assertEqual([ZERO, ZERO, HALF, HALF], chi(1,:), TOL, 'Chi' )
    @assertEqual([ZERO, ZERO, 4.0_defReal, ZERO], transOS(1,:), TOL, 'Transport XS O.S.' )
    @assertEqual([ZERO, ZERO, 5.5_defReal, ZERO], transFL(1,:), TOL, 'Transport XS F.L.' )
    @assertEqual([ZERO, ZERO, ZERO, ZERO, ZERO, TWO, ZERO, ZERO], P0(1,:), TOL, 'P0' )
    @assertEqual([ZERO, ZERO, ZERO, ZERO, ZERO, 1.5_defReal, ZERO, ZERO], P1(1,:), TOL, 'P1' )
    @assertEqual([ONE, ONE, ONE, ONE, ONE, TWO, ONE, ONE], prod(1,:), TOL, 'prod' )

    @assertEqual([ZERO, ZERO, ZERO, ZERO, ZERO, 0.6875_defReal, ZERO, ZERO], P2(1,:), TOL, 'P2' )
    @assertEqual([ZERO, ZERO, ZERO, ZERO, ZERO, -0.140625_defReal, ZERO, ZERO], P3(1,:), TOL, 'P3' )
    @assertEqual([ZERO, ZERO, ZERO, ZERO, ZERO, -0.7001953125_defReal, ZERO, ZERO], P4(1,:), TOL, 'P4' )
    @assertEqual([ZERO, ZERO, ZERO, ZERO, ZERO, -0.8327636719_defReal, ZERO, ZERO], P5(1,:), TOL, 'P5' )
    @assertEqual([ZERO, ZERO, ZERO, ZERO, ZERO, -0.5615539551_defReal, ZERO, ZERO], P6(1,:), TOL, 'P6' )
    @assertEqual([ZERO, ZERO, ZERO, ZERO, ZERO, -0.0683670044_defReal, ZERO, ZERO], P7(1,:), TOL, 'P7' )

    ! Test getting size
    @assertEqual(100, this % clerk_test1 % getSize(),'Test getSize():')

    ! Test correctness of output calls
    call out % init('dummyPrinter', fatalErrors = .false.)
    call this % clerk_test1 % print(out, mem)
    @assertTrue(out % isValid(), 'Test print():')

  end subroutine testScoring_clerk1

  !!
  !! Scoring test for clerk 2
  !!
@Test
  subroutine testScoring_clerk2(this)
    class(test_mgXsClerk), intent(inout) :: this
    character(:),allocatable             :: case
    type(scoreMemory)                    :: mem
    type(particle)                       :: p
    type(particleState)                  :: pFiss
    type(outputFile)                     :: out
    real(defReal), dimension(:,:), allocatable :: fiss, capt, transFL, transOS, &
                                                  nu, chi, P0, P1, prod
    real(defReal), parameter :: TOL = 1.0E-9

    ! Configure memory
    call mem % init(1000_longInt, 1)
    call this % clerk_test2 % setMemAddress(1_longInt)

    ! Report fission particles
    p % isMG = .false.
    p % w    = 0.5_defReal
    p % E    = 3.0_defReal
    call p % setMatIdx(2)

    ! Scoring
    pFiss = p
    call this % clerk_test2 % reportSpawn(p, pFiss, this % nucData, mem)

    p % E = 0.3_defReal

    ! Scoring
    pFiss = p
    call this % clerk_test2 % reportSpawn(p, pFiss, this % nucData, mem)

    ! Scoring
    call this % clerk_test2 % reportInColl(p, this % nucData, mem, .false.)

    p % preCollision % wgt = 0.2_defReal
    p % preCollision % E   = 0.3_defReal
    p % preCollision % isMG   = .false.
    p % E = 1.1_defReal

    call this % clerk_test2 % reportOutColl(p, N_2N, 0.75_defReal, this % nucData, mem)
    call mem % closeCycle(ONE)

    ! Process and get results
    call this % clerk_test2 % processRes(mem, capt, fiss, transFL, transOS, nu, chi, P0, P1, prod)

    ! Verify results of scoring
    @assertEqual([ZERO, ZERO, TWO], capt(1,:), TOL, 'Capture XS' )
    @assertEqual([ZERO, ZERO, 1.5_defReal], fiss(1,:), TOL, 'Fission XS' )
    @assertEqual([ZERO, ZERO, TWO], nu(1,:), TOL, 'NuFission XS' )
    @assertEqual([HALF, ZERO, HALF], chi(1,:), TOL, 'Chi' )
    @assertEqual([ZERO, ZERO, 4.0_defReal], transOS(1,:), TOL, 'Transport XS O.S.' )
    @assertEqual([ZERO, ZERO, 5.5_defReal], transFL(1,:), TOL, 'Transport XS F.L.' )
      @assertEqual([ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, TWO, ZERO], P0(1,:), TOL, 'P0' )
    @assertEqual([ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, 1.5_defReal, ZERO], P1(1,:), TOL, 'P1' )
    @assertEqual([ONE, ONE, ONE, ONE, ONE, ONE, ONE, TWO, ONE], prod(1,:), TOL, 'prod' )

    ! Test getting size
    @assertEqual(48, this % clerk_test2 % getSize(),'Test getSize():')

    ! Test correctness of output calls
    call out % init('dummyPrinter', fatalErrors = .false.)
    call this % clerk_test2 % print(out, mem)
    @assertTrue(out % isValid(), 'Test print():')

  end subroutine testScoring_clerk2

end module mgXsClerk_test
