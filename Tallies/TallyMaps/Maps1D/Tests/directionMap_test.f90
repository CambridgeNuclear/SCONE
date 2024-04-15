module directionMap_test
  use numPrecision
  use pFUnit_mod
  use particle_class,     only : particleState
  use dictionary_class,   only : dictionary
  use outputFile_class,   only : outputFile

  use directionMap_class, only : directionMap

  implicit none


@testCase
  type, extends(TestCase) :: test_directionMap
    private
    type(directionMap) :: map1
    type(directionMap) :: map2

  contains
    procedure :: setUp
    procedure :: tearDown

  end type test_directionMap

contains

  !!
  !! Sets up test_directionMap object
  !!
  subroutine setUp(this)
    class(test_directionMap), intent(inout) :: this
    type(dictionary)                        :: tempDict

    ! Build map with yz plane
    call tempDict % init(2)
    call tempDict % store('plane','yz')
    call tempDict % store('N', 4)

    call this % map1 % init(tempDict)
    call tempDict % kill()

    ! Build map with default plane
    call tempDict % init(1)
    call tempDict % store('N', 8)

    call this % map2 % init(tempDict)
    call tempDict % kill()

  end subroutine setUp

  !!
  !! Kills test_directionMap objects
  !!
  subroutine tearDown(this)
    class(test_directionMap), intent(inout) :: this

    call this % map1 % kill()
    call this % map2 % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test first map
  !!
@Test
  subroutine testMap1(this)
    class(test_directionMap), intent(inout)  :: this
    real(defReal),dimension(4),parameter     :: x   = [0.44_defReal, 15.8_defReal, 83.2_defReal, 999.1_defReal]
    real(defReal),dimension(4),parameter     :: phi = [-0.1_defReal, 0.84_defReal, 0.33_defReal, 0.64_defReal]*PI
    integer(shortInt),dimension(4),parameter :: RES_IDX = [2, 4, 3, 4]
    integer(shortInt),dimension(4)           :: idx
    type(particleState),dimension(4)         :: states

    ! Initialise states
    states(:) % dir(1) = x
    states(:) % dir(2) = cos(phi)
    states(:) % dir(3) = sin(phi)

    idx = this % map1 % map(states)

    @assertEqual(RES_IDX, idx)

  end subroutine testMap1

  !!
  !! Test second
  !!
@Test
  subroutine testMap2(this)
    class(test_directionMap), intent(inout)  :: this
    real(defReal),dimension(4),parameter     :: x = [0.9999_defReal, 0.87_defReal, -0.3_defReal, 0.18_defReal]
    real(defReal),dimension(4),parameter     :: y = [0.45_defReal, 0.88_defReal, 0.51_defReal, -0.92_defReal]
    real(defReal),dimension(4),parameter     :: z = [1.0_defReal, 39.8_defReal, 0.05_defReal, -12.2_defReal]
    integer(shortInt),dimension(4),parameter :: RES_IDX = [5, 6, 7, 3]
    integer(shortInt),dimension(4)           :: idx
    type(particleState),dimension(4)         :: states

    ! Initialise states
    states(:) % dir(1) = x
    states(:) % dir(2) = y
    states(:) % dir(3) = z

    idx = this % map2 % map(states)
    @assertEqual(RES_IDX, idx)

  end subroutine testMap2

  !!
  !! Test bin number retrival
  !!
@Test
  subroutine testBinNumber(this)
    class(test_directionMap), intent(inout) :: this

    ! Test that map is 1D
    @assertEqual(4, this % map1 % bins(0),'All bins')
    @assertEqual(4, this % map1 % bins(1),'1st dimension')
    @assertEqual(0, this % map2 % bins(2),'2nd dimension')
    @assertEqual(8, this % map2 % bins(1),'1st dimension')

    ! Get dimensionality
    @assertEqual(1, this % map1 % dimensions())
    @assertEqual(1, this % map2 % dimensions())

  end subroutine testBinNumber

  !!
  !! Test correctness of print subroutine
  !! Does not checks that values are correct, but that calls sequence is without errors
  !!
@Test
  subroutine testPrint(this)
    class(test_directionMap), intent(inout) :: this
    type(outputFile)                        :: out

    call out % init('dummyPrinter', fatalErrors = .false.)

    call this % map2 % print(out)
    @assertTrue(out % isValid(),'Radial map case')
    call out % reset()

    call this % map1 % print(out)
    @assertTrue(out % isValid(),'Unstruct map case')
    call out % reset()

  end subroutine testPrint


end module directionMap_test
