module sphericalMap_test
  use numPrecision
  use pFUnit_mod
  use particle_class,          only : particleState
  use dictionary_class,        only : dictionary
  use outputFile_class,        only : outputFile

  use sphericalMap_class,      only : sphericalMap

  implicit none


@testCase
  type, extends(TestCase) :: test_sphericalMap
    private
    type(sphericalMap) :: map_from_zero
    type(sphericalMap) :: map_from_min

  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_sphericalMap

contains

  !!
  !! Sets up test_sphericalMap object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_sphericalMap), intent(inout) :: this
    type(dictionary)                     :: tempDict

    ! Build map with default origin & minimum radius
    call tempDict % init(4)
    call tempDict % store('grid','lin')
    call tempDict % store('Rmax', 10.0_defReal)
    call tempDict % store('N', 20)

    call this % map_from_zero % init(tempDict)
    call tempDict % kill()

    ! Build map with diffrent origin & minimum radius
    call tempDict % init(5)
    call tempDict % store('origin', [ONE, ONE, ONE])
    call tempDict % store('grid', 'lin')
    call tempDict % store('Rmin', 5.0_defReal)
    call tempDict % store('Rmax', 10.0_defReal)
    call tempDict % store('N', 5)

    call this % map_from_min % init(tempDict)
    call tempDict % kill()

  end subroutine setUp

  !!
  !! Kills test_sphericalMap object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_sphericalMap), intent(inout) :: this

    call this % map_from_zero % kill()
    call this % map_from_min % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test default-initialised grid
  !!
@Test
  subroutine testFromOrigin(this)
    class(test_sphericalMap), intent(inout) :: this
    real(defReal),dimension(4),parameter     :: r = [0.4_defReal, 3.58_defReal, 8.9_defReal, 11.0_defReal]
    real(defReal), dimension(4),parameter    :: phi = [1.4_defReal, 3.98_defReal, 0.5_defReal, PI/2]
    real(defReal), dimension(4), parameter   :: tht = [ZERO, PI/2, PI/4, -PI/2]
    integer(shortInt),dimension(4),parameter :: RES_IDX = [1, 8, 18, 0]
    integer(shortInt),dimension(4)           :: idx
    type(particleState),dimension(4)         :: states

    ! Initialise states
    states(:) % r(1) = r * cos(phi) * sin(tht)
    states(:) % r(2) = r * sin(phi) * sin(tht)
    states(:) % r(3) = r * cos(tht)

    idx = this % map_from_zero % map(states)
    @assertEqual(RES_IDX, idx)

  end subroutine testFromOrigin

  !!
  !! Test grid with shifted origin & minimum radius
  !!
@Test
  subroutine testFromMin(this)
    class(test_sphericalMap), intent(inout) :: this
    real(defReal),dimension(4),parameter     :: r = [1.5_defReal, 5.5_defReal, 8.9_defReal, 11.0_defReal]
    real(defReal), dimension(4),parameter    :: phi = [1.4_defReal, 3.98_defReal, 0.5_defReal, PI/2]
    real(defReal), dimension(4), parameter   :: tht = [ZERO, PI/2, PI/4, -PI/2]
    integer(shortInt),dimension(4),parameter :: RES_IDX = [0, 1, 4, 0]
    integer(shortInt),dimension(4)           :: idx
    type(particleState),dimension(4)         :: states

    ! Initialise states
    states(:) % r(1) = r * cos(phi) * sin(tht)
    states(:) % r(2) = r * sin(phi) * sin(tht)
    states(:) % r(3) = r * cos(tht)

    ! Shift the origin
    states(:) % r(1) = states(:) % r(1) + ONE
    states(:) % r(2) = states(:) % r(2) + ONE
    states(:) % r(3) = states(:) % r(3) + ONE

    idx = this % map_from_min % map(states)
    @assertEqual(RES_IDX, idx)

  end subroutine testFromMin

  !!
  !! Test bin number retrival
  !!
@Test
  subroutine testBinNumber(this)
    class(test_sphericalMap), intent(inout) :: this

    @assertEqual(20, this % map_from_zero % bins(1),'1st Dimension')
    @assertEqual(20, this % map_from_zero % bins(0),'All bins')
    @assertEqual(0,  this % map_from_zero % bins(-3),'Invalid Dimension')

    ! Get dimensionality
    @assertEqual(1, this % map_from_zero % dimensions())

  end subroutine testBinNumber

  !!
  !! Test correctness of print subroutine
  !! Does not checks that values are correct, but that calls sequance is without errors
  !!
@Test
  subroutine testPrint(this)
    class(test_sphericalMap), intent(inout) :: this
    type(outputFile)                        :: out

    call out % init('dummyPrinter', fatalErrors = .false.)

    call this % map_from_zero % print(out)
    @assertTrue(out % isValid(),'Default-initialised map case')
    call out % reset()

    call this % map_from_min % print(out)
    @assertTrue(out % isValid(),'Map with minimum R')
    call out % reset()

  end subroutine testPrint


end module sphericalMap_test
