module endfTable_test
  use numPrecision
  use endfConstants
  use endfTable_class, only : endfTable
  use pFUnit_mod

  implicit none

contains

  !!
  !! Test table with lin-lin interpolation
  !!
@Test
  subroutine testLinLinTable()
    real(defReal),dimension(*),parameter :: x = [-1.0_defReal, 2.0_defReal, 2.5_defReal, 3.5_defReal]
    real(defReal),dimension(*),parameter :: y = [17.0_defReal, -2.0_defReal, 0.0_defReal, 1.5_defReal]
    type(endfTable)                      :: tab
    real(defReal),parameter :: TOL = 1.0E-6

    call tab % init(x, y)

    ! Check values at few points (including start and end of the table)
    @assertEqual(17.0_defReal, tab % at(-1.0_defReal), TOL)
    @assertEqual(1.5_defReal, tab % at(3.5_defReal), TOL)
    @assertEqual(-2.0_defReal, tab % at(2.0_defReal), TOL)

    ! Check interpolated values
    @assertEqual(32.0_defReal/3, tab % at(0.0_defReal), TOL)
    @assertEqual(0.75_defReal, tab % at(3.0_defReal), TOL)

  end subroutine testLinLinTable

  !!
  !! Test table with a single interpolation region that is not lin-lin
  !!
@Test
  subroutine testOneRegionInterpolation()
    real(defReal),dimension(*),parameter :: x = [-1.0_defReal, 2.0_defReal, 2.5_defReal, 3.5_defReal]
    real(defReal),dimension(*),parameter :: y = [17.0_defReal, -2.0_defReal, 0.0_defReal, 1.5_defReal]
    type(endfTable)                      :: tab
    real(defReal),parameter :: TOL = 1.0E-6

    ! Initialise with histogram interpolation
    call tab % init(x, y, [4], [histogramInterpolation])

    ! Check values at few points (including start and end of the table)
    @assertEqual(17.0_defReal, tab % at(-1.0_defReal), TOL)
    @assertEqual(0.0_defReal, tab % at(3.5_defReal), TOL)
    @assertEqual(-2.0_defReal, tab % at(2.0_defReal), TOL)

    ! Check interpolated values
    @assertEqual(17.0_defReal, tab % at(0.0_defReal), TOL)
    @assertEqual(0.0_defReal, tab % at(3.0_defReal), TOL)

  end subroutine testOneRegionInterpolation

  !!
  !! Test table with multiple interpolation regions
  !!
@Test
  subroutine testMultipleInterRegions()
    real(defReal),dimension(*),parameter :: x = [-1.0_defReal, 2.0_defReal, 2.5_defReal, 2.5_defReal, 3.5_defReal]
    real(defReal),dimension(*),parameter :: y = [17.0_defReal, 2.0_defReal, 0.0_defReal, 1.0_defReal, 1.5_defReal]
    type(endfTable)                      :: tab
    real(defReal),parameter :: TOL = 1.0E-6

    ! Initialise with linLog, logLin and logLog interpolations
    call tab % init(x, y, [2,4,5], [loglinInterpolation, linLogInterpolation, logLogInterpolation])

    ! Check values at few points (including start and end of the table)
    @assertEqual(17.0_defReal, tab % at(-1.0_defReal), TOL)
    @assertEqual(1.5_defReal, tab % at(3.5_defReal), TOL)
    @assertEqual(2.0_defReal, tab % at(2.0_defReal), TOL)

    ! Check y-log x-lin interpolation
    @assertEqual(8.329954_defReal ,tab % at(0.0_defReal), TOL)

    ! Check y-lin x-log interpolation
    @assertEqual(1.5627016_defReal, tab % at(2.1_defReal),TOL)

    ! Check y-log x-log interpolation
    @assertEqual(1.346458993_defReal, tab % at(3.2_defReal), TOL)

  end subroutine testMultipleInterRegions
    
end module endfTable_test
