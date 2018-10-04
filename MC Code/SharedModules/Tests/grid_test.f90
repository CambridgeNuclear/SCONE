module grid_test
  use numPrecision
  use universalVariables
  use pfUnit_mod
  use grid_class, only : grid

  implicit none

@TestCase
  type, extends(TestCase) :: test_grid
    type(grid) :: linGrid
    type(grid) :: logGrid
    type(grid) :: unstructGrid
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_grid

  !! Patameters
  real(defReal), parameter :: FP_TOL = 50*epsilon(1.0_defReal)

contains

  !!
  !! Initialise test case for structured linear grid
  !!
  subroutine setUp(this)
    class(test_grid), intent(inout) :: this
    character(nameLen)              :: type

    ! Typical PWR Assembly Grid - Linear Grid
    type = 'lin'
    call this % linGrid % init(-10.71_defReal, 10.71_defReal, 17, type)

    ! Equilethargic 70 Group Spectrum - Log Grid
    type = 'log'
    call this % logGrid % init(1.0E-11_defReal, 20.0_defReal, 70, type)

    ! CASMO 8 Group Structure - Unstructured grid
    call this % unstructGrid % init([1.00e-11_defReal, &
                                     5.80e-08_defReal, &
                                     1.40e-07_defReal, &
                                     2.80e-07_defReal, &
                                     6.25e-07_defReal, &
                                     4.00e-06_defReal, &
                                     0.005530_defReal, &
                                     0.821000_defReal, &
                                     10.00000_defReal])

  end subroutine setUp

  !!
  !! Deconstruct test case for structured linear grid
  !!
  subroutine tearDown(this)
    class(test_grid), intent(inout) :: this

    ! Do nothing *** TODO: Add kill later

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!!Tests on linear structured grid
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

!!
!! Test valid quary about the bin boundary
!!
@test
  subroutine testValidBinQuery_lin(this)
    class(test_grid), intent(inout)          :: this
    real(defReal),dimension(5)               :: bin
    integer(shortInt),dimension(5),parameter :: idx = [1,4,7,14,18]

    bin = this % linGrid % bin(idx)
    @assertEqual(-10.71_defReal, bin(1), abs(FP_TOL*bin(1)) )
    @assertEqual(-6.93_defReal,  bin(2), abs(FP_TOL*bin(2)) )
    @assertEqual(-3.15_defReal,  bin(3), abs(FP_TOL*bin(3)) )
    @assertEqual(5.67_defReal,   bin(4), abs(FP_TOL*bin(4)) )
    @assertEqual(10.71_defReal,  bin(5), abs(FP_TOL*bin(5)) )

  end subroutine testValidBinQuery_lin

!!
!! Test invalid quary about the bin boundary
!!
@test
  subroutine testInvalidBinQuery_lin(this)
    class(test_grid), intent(inout)          :: this
    real(defReal),dimension(3)               :: bin
    integer(shortInt),dimension(3),parameter :: idx = [-1,0, 20]
    real(defReal),parameter                  :: INVALID_BIN = -huge(bin)

    bin = this % linGrid % bin(idx)
    @assertEqual(INVALID_BIN, bin(1))
    @assertEqual(INVALID_BIN, bin(2))
    @assertEqual(INVALID_BIN, bin(3))

  end subroutine testInvalidBinQuery_lin

!!
!! Test search
!!
@test
  subroutine testSearch_lin(this)
    class(test_grid), intent(inout) :: this
    integer(shortInt),dimension(5)  :: idx
    real(defReal),dimension(5)      :: keys

    keys = [-10.71_defReal, 3.13245_defReal, -8.96_defReal, -20.0_defReal, 10.72_defReal]
    idx = this % linGrid % search(keys)

    @assertEqual(1,  idx(1))
    @assertEqual(11, idx(2))
    @assertEqual(2,  idx(3))

    ! Invalid Searches
    @assertEqual(valueOutsideArray, idx(4))
    @assertEqual(valueOutsideArray, idx(5))

  end subroutine testSearch_lin

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!!Tests on logarithmic structured grid
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

!!
!! Test valid quary about the bin boundary
!!
@test
  subroutine testValidBinQuery_log(this)
    class(test_grid), intent(inout)          :: this
    real(defReal),dimension(4)               :: bin
    integer(shortInt),dimension(4),parameter :: idx = [1,35,70,71]

    bin = this % logGrid % bin(idx)
    @assertEqual(1.0E-11_defReal,               bin(1), abs(FP_TOL*bin(1)) )
    @assertEqual(9.435957972757912E-6_defReal,  bin(2), abs(FP_TOL*bin(2)) )
    @assertEqual(13.344459739056777_defReal,    bin(3), abs(FP_TOL*bin(3)) )
    @assertEqual(20.0_defReal,                  bin(4), abs(FP_TOL*bin(4)) )


  end subroutine testValidBinQuery_log

!!
!! Test invalid quary about the bin boundary
!!
@test
  subroutine testInvalidBinQuery_log(this)
    class(test_grid), intent(inout)          :: this
    real(defReal),dimension(3)               :: bin
    integer(shortInt),dimension(3),parameter :: idx = [-1,0, 72]
    real(defReal),parameter                  :: INVALID_BIN = -huge(bin)

    bin = this % logGrid % bin(idx)
    @assertEqual(INVALID_BIN, bin(1))
    @assertEqual(INVALID_BIN, bin(2))
    @assertEqual(INVALID_BIN, bin(3))

  end subroutine testInvalidBinQuery_log

!!
!! Test search
!!
@test
  subroutine testSearch_log(this)
    class(test_grid), intent(inout) :: this
    integer(shortInt),dimension(5)  :: idx
    real(defReal),dimension(5)      :: keys

    keys = [1.0E-11_defReal, 6.7E-4_defReal, 1.0_defReal, 9.0E-12_defReal, 21.0_defReal]
    idx = this % logGrid % search(keys)

    @assertEqual(1,  idx(1))
    @assertEqual(45, idx(2))
    @assertEqual(63, idx(3))

    ! Invalid Searches
    @assertEqual(valueOutsideArray, idx(4))
    @assertEqual(valueOutsideArray, idx(5))

  end subroutine testSearch_log

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!!Tests on unstructured grid
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

!!
!! Test valid quary about the bin boundary
!!
@test
  subroutine testValidBinQuery_unstruct(this)
    class(test_grid), intent(inout)          :: this
    real(defReal),dimension(4)               :: bin
    integer(shortInt),dimension(4),parameter :: idx = [1,3,7,9]

    bin = this % unstructGrid % bin(idx)
    @assertEqual(1.0E-11_defReal,  bin(1), abs(FP_TOL*bin(1)) )
    @assertEqual(1.4E-7_defReal,   bin(2), abs(FP_TOL*bin(2)) )
    @assertEqual(0.00553_defReal,  bin(3), abs(FP_TOL*bin(3)) )
    @assertEqual(10.0_defReal,     bin(4), abs(FP_TOL*bin(4)) )

  end subroutine testValidBinQuery_unstruct

!!
!! Test invalid quary about the bin boundary
!!
@test
  subroutine testInvalidBinQuery_unstruct(this)
    class(test_grid), intent(inout)          :: this
    real(defReal),dimension(3)               :: bin
    integer(shortInt),dimension(3),parameter :: idx = [-1,0, 72]
    real(defReal),parameter                  :: INVALID_BIN = -huge(bin)

    bin = this % unstructGrid % bin(idx)
    @assertEqual(INVALID_BIN, bin(1))
    @assertEqual(INVALID_BIN, bin(2))
    @assertEqual(INVALID_BIN, bin(3))

  end subroutine testInvalidBinQuery_unstruct

!!
!! Test search
!!
@test
  subroutine testSearch_unstruct(this)
    class(test_grid), intent(inout) :: this
    integer(shortInt),dimension(5)  :: idx
    real(defReal),dimension(5)      :: keys

    keys = [1.0E-11_defReal, 6.7E-4_defReal, 1.0_defReal, 9.0E-12_defReal, 21.0_defReal]
    idx = this % unstructGrid % search(keys)

    @assertEqual(1,  idx(1))
    @assertEqual(6,  idx(2))
    @assertEqual(8,  idx(3))

    ! Invalid Searches
    @assertEqual(valueOutsideArray, idx(4))
    @assertEqual(valueOutsideArray, idx(5))

  end subroutine testSearch_unstruct

end module grid_test
