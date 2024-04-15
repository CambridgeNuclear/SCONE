module cylindricalRadMap_test
  use numPrecision
  use pFUnit_mod
  use particle_class,          only : particleState
  use dictionary_class,        only : dictionary
  use outputFile_class,        only : outputFile

  use cylindricalRadMap_class, only : cylindricalRadMap

  implicit none


@testCase
  type, extends(TestCase) :: test_cylindricalRadMap
    private
    type(cylindricalRadMap) :: map_radial
    type(cylindricalRadMap) :: map_unstruct
    type(cylindricalRadMap) :: map_equivol

  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_cylindricalRadMap

contains

  !!
  !! Sets up test_cylindricalRadMap object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_cylindricalRadMap), intent(inout) :: this
    type(dictionary)                             :: tempDict

    ! Build radial map with different orientation & minimum radius
    call tempDict % init(5)
    call tempDict % store('orientation','x')
    call tempDict % store('rGrid','equivolume')
    call tempDict % store('Rmin', 2.0_defReal)
    call tempDict % store('Rmax', 10.0_defReal)
    call tempDict % store('rN', 5)

    call this % map_radial % init(tempDict)
    call tempDict % kill()

    ! Build map with different origin & unstruct bins
    call tempDict % init(3)
    call tempDict % store('origin',[ONE, ONE])
    call tempDict % store('rGrid','unstruct')
    call tempDict % store('bins', [1.5_defReal, 2.3_defReal, 3.8_defReal, 8.0_defReal])

    call this % map_unstruct % init(tempDict)
    call tempDict % kill()

    ! Build map with equivolume bins
    call tempDict % init(3)
    call tempDict % store('rGrid','equivolume')
    call tempDict % store('Rmax', 8.0_defReal)
    call tempDict % store('rN', 4)

    call this % map_equivol % init(tempDict)
    call tempDict % kill()

  end subroutine setUp

  !!
  !! Kills test_cylindricalRadMap object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_cylindricalRadMap), intent(inout) :: this

    call this % map_radial % kill()
    call this % map_unstruct % kill()
    call this % map_equivol % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test radial map with different orientation & minimum radius
  !!
@Test
  subroutine testRadial(this)
    class(test_cylindricalRadMap), intent(inout) :: this
    real(defReal),dimension(4),parameter         :: r = [0.4_defReal, 5.38_defReal, 8.9_defReal, 9.1_defReal]
    real(defReal), dimension(4),parameter        :: phi = [1.4_defReal, 3.0_defReal, 0.5_defReal, PI/2]
    real(defReal), dimension(4),parameter        :: z = [1.0_defReal, 39.8_defReal, 0.05_defReal, -12.2_defReal]
    integer(shortInt),dimension(4),parameter     :: RES_IDX = [0, 2, 4, 5]
    integer(shortInt),dimension(4)               :: idx
    type(particleState),dimension(4)             :: states

    ! Initialise states
    states(:) % r(1) = z
    states(:) % r(2) = r * cos(phi)
    states(:) % r(3) = r * sin(phi)

    idx = this % map_radial % map(states)

    @assertEqual(RES_IDX, idx)

  end subroutine testRadial

  !!
  !! Test map with different origin & unstruct bins
  !!
@Test
  subroutine testUnstruct(this)
    class(test_cylindricalRadMap), intent(inout) :: this
    real(defReal),dimension(4),parameter         :: r = [1.52_defReal, 5.5_defReal, 8.9_defReal, 2.88_defReal]
    real(defReal), dimension(4),parameter        :: phi = [1.4_defReal, 3.0_defReal, 0.5_defReal, PI/2]
    real(defReal), dimension(4), parameter       :: z = [1.0_defReal, 39.8_defReal, 0.05_defReal, -12.2_defReal]
    integer(shortInt),dimension(4),parameter     :: RES_IDX = [1, 3, 0, 2]
    integer(shortInt),dimension(4)               :: idx
    type(particleState),dimension(4)             :: states

    ! Initialise states
    states(:) % r(1) = r * cos(phi)
    states(:) % r(2) = r * sin(phi)
    states(:) % r(3) = z

    ! Shift the origin
    states(:) % r(1) = states(:) % r(1) + ONE
    states(:) % r(2) = states(:) % r(2) + ONE
    states(:) % r(3) = states(:) % r(3)

    idx = this % map_unstruct % map(states)
    @assertEqual(RES_IDX, idx)

  end subroutine testUnstruct

  !!
  !! Test map with equivolume bins
  !!
@Test
  subroutine testEquivol(this)
    class(test_cylindricalRadMap), intent(inout) :: this
    real(defReal),dimension(4),parameter         :: r = [1.82_defReal, 4.68_defReal, 7.9_defReal, 17.01_defReal]
    real(defReal), dimension(4),parameter        :: phi = [1.4_defReal, 3.0_defReal, 0.5_defReal, PI/2]
    real(defReal), dimension(4), parameter       :: z = [1.0_defReal, 39.8_defReal, 0.05_defReal, -12.2_defReal]
    integer(shortInt),dimension(4),parameter     :: RES_IDX = [1, 2, 4, 0]
    integer(shortInt),dimension(4)               :: idx
    type(particleState),dimension(4)             :: states

    ! Initialise states
    states(:) % r(1) = r * cos(phi)
    states(:) % r(2) = r * sin(phi)
    states(:) % r(3) = z

    idx = this % map_equivol % map(states)
    @assertEqual(RES_IDX, idx)

  end subroutine testEquivol

  !!
  !! Test bin number retrival
  !!
@Test
  subroutine testBinNumber(this)
    class(test_cylindricalRadMap), intent(inout) :: this

    ! Test that map is 1D
    @assertEqual(5, this % map_radial % bins(0),'All bins')
    @assertEqual(5, this % map_radial % bins(1),'1st dimension')
    @assertEqual(0, this % map_radial % bins(2),'2nd dimension')
    @assertEqual(3, this % map_unstruct % bins(1),'1st dimension')
    @assertEqual(4, this % map_equivol % bins(1),'1st dimension')

    ! Get dimensionality
    @assertEqual(1, this % map_radial % dimensions())
    @assertEqual(1, this % map_unstruct % dimensions())

  end subroutine testBinNumber

  !!
  !! Test correctness of print subroutine
  !! Does not checks that values are correct, but that calls sequence is without errors
  !!
@Test
  subroutine testPrint(this)
    class(test_cylindricalRadMap), intent(inout) :: this
    type(outputFile)                             :: out

    call out % init('dummyPrinter', fatalErrors = .false.)

    call this % map_radial % print(out)
    @assertTrue(out % isValid(),'Radial map case')
    call out % reset()

    call this % map_unstruct % print(out)
    @assertTrue(out % isValid(),'Unstruct map case')
    call out % reset()

    call this % map_equivol % print(out)
    @assertTrue(out % isValid(),'Unstruct map case')
    call out % reset()

  end subroutine testPrint


end module cylindricalRadMap_test
