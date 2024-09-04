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
    call tempDict % init(3)
    call tempDict % store('N', 9)
    call tempDict % store('min', 90)
    call tempDict % store('max', 180)

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
    real(defReal),dimension(4),parameter     :: phi = [100.1_defReal, 20.84_defReal, 333.9_defReal, 264.8_defReal]*PI/180.0_defReal
    integer(shortInt),dimension(4),parameter :: RES_IDX = [4, 3, 2, 1]
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
    real(defReal),dimension(4),parameter     :: z   = [0.44_defReal, 15.8_defReal, 83.2_defReal, 999.1_defReal]
    real(defReal),dimension(4),parameter     :: phi = [170.1_defReal, 90.84_defReal, 133.9_defReal, 264.8_defReal]*PI/180.0_defReal
    integer(shortInt),dimension(4),parameter :: RES_IDX = [9, 1, 5, 0]
    integer(shortInt),dimension(4)           :: idx
    type(particleState),dimension(4)         :: states

    ! Initialise states
    states(:) % dir(1) = cos(phi)
    states(:) % dir(2) = sin(phi)
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
    @assertEqual(9, this % map2 % bins(1),'1st dimension')

    ! Get dimensionality
    @assertEqual(1, this % map1 % dimensions())
    @assertEqual(1, this % map2 % dimensions())

  end subroutine testBinNumber

  !!
  !! Test correctness of print subroutine
  !! Does not check that values are correct, but that call sequence is without errors
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
