module energyGrid_test
  use numPrecision
  use universalVariables
  use pfUnit_mod
  use energyGrid_class, only : energyGrid

  implicit none

@TestCase
  type, extends(TestCase) :: test_energyGrid
    type(energyGrid) :: linGrid
    type(energyGrid) :: logGrid
    type(energyGrid) :: unstructGrid
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_energyGrid

  !! Patameters
  real(defReal), parameter :: FP_TOL = 50*epsilon(1.0_defReal)

contains

  !!
  !! Initialise test case for structured linear energyGrid
  !!
  subroutine setUp(this)
    class(test_energyGrid), intent(inout) :: this
    character(nameLen)                    :: type

    ! Simple linear energy grid
    type = 'lin'
    call this % linGrid % init(ZERO, 20.0_defReal, 20, type)

    ! Simple logarithmic (equilethargic) grid
    type = 'log'
    call this % logGrid % init(1.0E-9_defReal, 1.0_defReal, 9, type)

    ! Simple unstructured grid CASMO-8
    call this % unstructGrid % init([10.00000_defReal, &
                                     0.821000_defReal, &
                                     0.005530_defReal, &
                                     4.00e-06_defReal, &
                                     6.25e-07_defReal, &
                                     2.80e-07_defReal, &
                                     1.40e-07_defReal, &
                                     5.80e-08_defReal, &
                                     1.00e-11_defReal])

  end subroutine setUp

  !!
  !! Deconstruct test case for structured linear energyGrid
  !!
  subroutine tearDown(this)
    class(test_energyGrid), intent(inout) :: this

    ! Do nothing *** TODO: Add kill later

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Proper Tests begin here
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test searching linear grid
  !!
@Test
  subroutine testSearchLinEnergyGrid(this)
    class(test_energyGrid), intent(inout) :: this

    ! Valid searches
    @assertEqual(1, this % linGrid % search(19.5_defReal))
    @assertEqual(20, this % linGrid % search(0.5_defReal))
    @assertEqual(11, this % linGrid % search(9.9_defReal))

    ! Invalid searches
    @assertEqual(valueOutsideArray, this % linGrid % search(20.5_defReal))
    @assertEqual(valueOutsideArray, this % linGrid % search(-0.1_defReal))
    @assertEqual(valueOutsideArray, this % linGrid % search(-2.1_defReal))

  end subroutine testSearchLinEnergyGrid

  !!
  !! Test searching logarithmic grid
  !!
@Test
  subroutine testSearchLogEnergyGrid(this)
    class(test_energyGrid), intent(inout) :: this

    ! Valid searches
    @assertEqual(1, this % logGrid % search(0.7_defReal))
    @assertEqual(6, this % logGrid % search(1.1E-6_defReal))
    @assertEqual(9, this % logGrid % search(2.0E-9_defReal))

    ! Invalid searches
    @assertEqual(valueOutsideArray, this % logGrid % search(1.1_defReal))
    @assertEqual(valueOutsideArray, this % logGrid % search(-0.1_defReal))
    @assertEqual(valueOutsideArray, this % logGrid % search(1.0E-11_defReal))

  end subroutine testSearchLogEnergyGrid

  !!
  !! Test searching unstructured grid
  !!
@Test
  subroutine testSearchUnstructEnergyGrid(this)
    class(test_energyGrid), intent(inout) :: this

    ! Valid searches
    @assertEqual(1, this % unstructGrid % search(1.7_defReal))
    @assertEqual(3, this % unstructGrid % search(5.1E-5_defReal))
    @assertEqual(8, this % unstructGrid % search(2.0E-9_defReal))

    ! Invalid searches
    @assertEqual(valueOutsideArray, this % unstructGrid % search(11.0_defReal))
    @assertEqual(valueOutsideArray, this % unstructGrid % search(-0.1_defReal))
    @assertEqual(valueOutsideArray, this % unstructGrid % search(0.9E-11_defReal))

  end subroutine testSearchUnstructEnergyGrid

  !!
  !! Test getting bin values from linear grid
  !!
@Test
  subroutine testBinLin(this)
    class(test_energyGrid), intent(inout) :: this
    real(defReal),parameter :: TOL =1.0E-9

    ! Valid Bins
    @assertEqual(20.0_defReal, this % linGrid % bin(1), TOL * 20.0_defReal)
    @assertEqual(0.0_defReal, this % linGrid % bin(21), TOL)
    @assertEqual(17.0_defReal, this % linGrid % bin(4), TOL * 17.0_defReal)

    ! Invalid Bins
    @assertEqual(-huge(ONE), this % linGrid % bin(22))
    @assertEqual(-huge(ONE), this % linGrid % bin(0))

  end subroutine testBinLin

  !!
  !! Test getting bin values from logarithmic grid
  !!
@Test
  subroutine testBinLog(this)
    class(test_energyGrid), intent(inout) :: this
    real(defReal),parameter :: TOL =1.0E-9

    ! Valid bins
    @assertEqual(1.0_defReal, this % logGrid % bin(1), TOL * 1.0_defReal)
    @assertEqual(1.0E-8_defReal, this % logGrid % bin(9), TOL * 1.0E-8_defReal)
    @assertEqual(1.0E-9_defReal, this % logGrid % bin(10), TOL * 1.0E-9_defReal)
    @assertEqual(1.0E-3_defReal, this % logGrid % bin(4), TOL * 1.0E-3_defReal)

    ! Invalid Bins
    @assertEqual(-huge(ONE), this % logGrid % bin(11))
    @assertEqual(-huge(ONE), this % logGrid % bin(0))

  end subroutine testBinLog

  !!
  !! Test getting bin values from logarithmic grid
  !!
@Test
  subroutine testBinUnstruct(this)
    class(test_energyGrid), intent(inout) :: this
    real(defReal),parameter :: TOL =1.0E-9

    ! Valid bins
    @assertEqual(10.00000_defReal, this % unstructGrid % bin(1), TOL * 10.00000_defReal)
    @assertEqual(5.80e-08_defReal, this % unstructGrid % bin(8), TOL * 5.80e-08_defReal)
    @assertEqual(1.00e-11_defReal, this % unstructGrid % bin(9), TOL * 1.00e-11_defReal)
    @assertEqual(4.00e-06_defReal, this % unstructGrid % bin(4), TOL * 4.00e-06_defReal)

    ! Invalid Bins
    @assertEqual(-huge(ONE), this % unstructGrid % bin(10))
    @assertEqual(-huge(ONE), this % unstructGrid % bin(0))

  end subroutine testBinUnstruct

  !!
  !! Test getting size (number of bins/groups)
  !!
@Test
  subroutine testGetSize(this)
    class(test_energyGrid), intent(inout) :: this

    @assertEqual(20, this % linGrid % getSize())
    @assertEqual(9, this % logGrid % getSize())
    @assertEqual(8, this % unstructGrid % getSize())

  end subroutine testGetSize

end module energyGrid_test
