module radialMap_test

  use numPrecision
  use funit
  use particle_class,   only : particleState
  use dictionary_class, only : dictionary
  use outputFile_class, only : outputFile

  use radialMap_class,  only : radialMap

  implicit none


@testCase
  type, extends(TestCase) :: test_radialMap
    private
    type(radialMap) :: map_cyl_linear
    type(radialMap) :: map_cyl_equivol
    type(radialMap) :: map_cyl_unstruct
    type(radialMap) :: map_sph_from_zero
    type(radialMap) :: map_sph_from_min
    type(radialMap) :: map_sph_equivol

  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_radialMap

contains

  !!
  !! Sets up test_radialMap object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_radialMap), intent(inout) :: this
    type(dictionary)                     :: tempDict

    ! Build cylindrical map with linear bins
    call tempDict % init(4)
    call tempDict % store('axis','x')
    call tempDict % store('grid','lin')
    call tempDict % store('max', 8.0_defReal)
    call tempDict % store('N', 4)

    call this % map_cyl_linear % init(tempDict)
    call tempDict % kill()

    ! Build cylindrical map with different orientation & minimum radius
    call tempDict % init(5)
    call tempDict % store('axis','z')
    call tempDict % store('grid','equivolume')
    call tempDict % store('min', 2.0_defReal)
    call tempDict % store('max', 10.0_defReal)
    call tempDict % store('N', 5)

    call this % map_cyl_equivol % init(tempDict)
    call tempDict % kill()

    ! Build cylindrical map with different origin & unstruct bins
    call tempDict % init(4)
    call tempDict % store('axis','z')
    call tempDict % store('origin',[ONE, ONE, ZERO])
    call tempDict % store('grid','unstruct')
    call tempDict % store('bins', [1.5_defReal, 2.3_defReal, 3.8_defReal, 8.0_defReal])

    call this % map_cyl_unstruct % init(tempDict)
    call tempDict % kill()

    ! Build spherical map with default origin & minimum radius
    call tempDict % init(3)
    call tempDict % store('grid','lin')
    call tempDict % store('max', 10.0_defReal)
    call tempDict % store('N', 20)

    call this % map_sph_from_zero % init(tempDict)
    call tempDict % kill()

    ! Build spherical map with diffrent origin & minimum radius
    call tempDict % init(5)
    call tempDict % store('origin', [ONE, ONE, ONE])
    call tempDict % store('grid', 'lin')
    call tempDict % store('min', 5.0_defReal)
    call tempDict % store('max', 10.0_defReal)
    call tempDict % store('N', 5)

    call this % map_sph_from_min % init(tempDict)
    call tempDict % kill()

    ! Build spherical map with equivolume bins
    call tempDict % init(4)
    call tempDict % store('grid', 'equivolume')
    call tempDict % store('min', 2.0_defReal)
    call tempDict % store('max', 20.0_defReal)
    call tempDict % store('N', 8)

    call this % map_sph_equivol % init(tempDict)
    call tempDict % kill()

  end subroutine setUp

  !!
  !! Kills test_radialMap object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_radialMap), intent(inout) :: this

    call this % map_cyl_linear % kill()
    call this % map_cyl_equivol % kill()
    call this % map_cyl_unstruct % kill()
    call this % map_sph_from_zero % kill()
    call this % map_sph_from_min % kill()
    call this % map_sph_equivol % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test cylindrical map with different orientation & minimum radius
  !!
@Test
  subroutine testCylLinear(this)
    class(test_radialMap), intent(inout)     :: this
    real(defReal),dimension(4),parameter     :: r = [0.4_defReal, 5.38_defReal, 7.9_defReal, 9.1_defReal]
    real(defReal), dimension(4),parameter    :: phi = [1.4_defReal, 3.0_defReal, 0.5_defReal, PI/2]
    real(defReal), dimension(4),parameter    :: z = [1.0_defReal, 39.8_defReal, 0.05_defReal, -12.2_defReal]
    integer(shortInt),dimension(4),parameter :: RES_IDX = [1, 3, 4, 0]
    integer(shortInt),dimension(4)           :: idx
    type(particleState),dimension(4)         :: states

    ! Initialise states
    states(:) % r(1) = z
    states(:) % r(2) = r * cos(phi)
    states(:) % r(3) = r * sin(phi)

    idx = this % map_cyl_linear % map(states)

    @assertEqual(RES_IDX, idx)

  end subroutine testCylLinear

  !!
  !! Test cylindrical map with equivolume bins
  !!
@Test
  subroutine testCylEquivol(this)
    class(test_radialMap), intent(inout)     :: this
    real(defReal),dimension(4),parameter     :: r = [1.82_defReal, 4.68_defReal, 7.9_defReal, 9.01_defReal]
    real(defReal), dimension(4),parameter    :: phi = [1.4_defReal, 3.0_defReal, 0.5_defReal, PI/2]
    real(defReal), dimension(4), parameter   :: z = [1.0_defReal, 39.8_defReal, 0.05_defReal, -12.2_defReal]
    integer(shortInt),dimension(4),parameter :: RES_IDX = [0, 1, 4, 5]
    integer(shortInt),dimension(4)           :: idx
    type(particleState),dimension(4)         :: states

    ! Initialise states
    states(:) % r(1) = r * cos(phi)
    states(:) % r(2) = r * sin(phi)
    states(:) % r(3) = z

    idx = this % map_cyl_equivol % map(states)
    @assertEqual(RES_IDX, idx)

  end subroutine testCylEquivol

  !!
  !! Test cylindrical map with different origin & unstruct bins
  !!
@Test
  subroutine testCylUnstruct(this)
    class(test_radialMap), intent(inout)     :: this
    real(defReal),dimension(4),parameter     :: r = [1.52_defReal, 5.5_defReal, 8.9_defReal, 2.88_defReal]
    real(defReal), dimension(4),parameter    :: phi = [1.4_defReal, 3.0_defReal, 0.5_defReal, PI/2]
    real(defReal), dimension(4), parameter   :: z = [1.0_defReal, 39.8_defReal, 0.05_defReal, -12.2_defReal]
    integer(shortInt),dimension(4),parameter :: RES_IDX = [1, 3, 0, 2]
    integer(shortInt),dimension(4)           :: idx
    type(particleState),dimension(4)         :: states

    ! Initialise states
    states(:) % r(1) = r * cos(phi)
    states(:) % r(2) = r * sin(phi)
    states(:) % r(3) = z

    ! Shift the origin
    states(:) % r(1) = states(:) % r(1) + ONE
    states(:) % r(2) = states(:) % r(2) + ONE
    states(:) % r(3) = states(:) % r(3)

    idx = this % map_cyl_unstruct % map(states)
    @assertEqual(RES_IDX, idx)

  end subroutine testCylUnstruct

  !!
  !! Test spherical map with default-initialised grid
  !!
@Test
  subroutine testSphFromOrigin(this)
    class(test_radialMap), intent(inout)     :: this
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

    idx = this % map_sph_from_zero % map(states)
    @assertEqual(RES_IDX, idx)

  end subroutine testSphFromOrigin

  !!
  !! Test spherical map with grid with shifted origin & minimum radius
  !!
@Test
  subroutine testSphFromMin(this)
    class(test_radialMap), intent(inout)     :: this
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

    idx = this % map_sph_from_min % map(states)
    @assertEqual(RES_IDX, idx)

  end subroutine testSphFromMin

  !!
  !! Test spherical map with grid with equivolume bins
  !!
@Test
  subroutine testSphEquivol(this)
    class(test_radialMap), intent(inout)     :: this
    real(defReal),dimension(4),parameter     :: r = [1.5_defReal, 5.5_defReal, 18.9_defReal, 11.0_defReal]
    real(defReal), dimension(4),parameter    :: phi = [1.4_defReal, 3.98_defReal, 0.5_defReal, PI/2]
    real(defReal), dimension(4), parameter   :: tht = [ZERO, PI/2, PI/4, -PI/2]
    integer(shortInt),dimension(4),parameter :: RES_IDX = [0, 1, 7, 2]
    integer(shortInt),dimension(4)           :: idx
    type(particleState),dimension(4)         :: states

    ! Initialise states
    states(:) % r(1) = r * cos(phi) * sin(tht)
    states(:) % r(2) = r * sin(phi) * sin(tht)
    states(:) % r(3) = r * cos(tht)

    idx = this % map_sph_equivol % map(states)
    @assertEqual(RES_IDX, idx)

  end subroutine testSphEquivol

  !!
  !! Test bin number retrival
  !!
@Test
  subroutine testBinNumber(this)
    class(test_radialMap), intent(inout) :: this

    ! Test that map is 1D
    @assertEqual(4, this % map_cyl_linear % bins(0),'All bins')
    @assertEqual(4, this % map_cyl_linear % bins(1),'1st dimension')
    @assertEqual(0, this % map_cyl_linear % bins(2),'2nd dimension')
    @assertEqual(3, this % map_cyl_unstruct % bins(1),'1st dimension')
    @assertEqual(5, this % map_cyl_equivol % bins(1),'1st dimension')
    @assertEqual(20, this % map_sph_from_zero % bins(1),'1st Dimension')
    @assertEqual(20, this % map_sph_from_zero % bins(0),'All bins')
    @assertEqual(0,  this % map_sph_from_min % bins(2),'Invalid Dimension')

    ! Get dimensionality
    @assertEqual(1, this % map_cyl_linear % dimensions())
    @assertEqual(1, this % map_cyl_unstruct % dimensions())
    @assertEqual(1, this % map_sph_from_min % dimensions())

  end subroutine testBinNumber

  !!
  !! Test correctness of print subroutine
  !! Does not check that values are correct, but that call sequence is without errors
  !!
@Test
  subroutine testPrint(this)
    class(test_radialMap), intent(inout) :: this
    type(outputFile)                     :: out

    call out % init('dummyPrinter', fatalErrors = .false.)

    call this % map_cyl_linear % print(out)
    @assertTrue(out % isValid(),'Linear map case (cylindrical)')
    call out % reset()

    call this % map_cyl_equivol % print(out)
    @assertTrue(out % isValid(),'Equivolume map case (cylindrical)')
    call out % reset()

    call this % map_cyl_unstruct % print(out)
    @assertTrue(out % isValid(),'Unstruct map case (cylindrical)')
    call out % reset()

    call this % map_sph_from_zero % print(out)
    @assertTrue(out % isValid(),'Linear map case from zero (spherical)')
    call out % reset()

    call this % map_sph_from_min % print(out)
    @assertTrue(out % isValid(),'Linear map case from minimum radius (spherical)')
    call out % reset()

    call this % map_sph_equivol % print(out)
    @assertTrue(out % isValid(),'Equivolume map case (spherical)')
    call out % reset()

  end subroutine testPrint


end module radialMap_test
