module endfTable_test
  use numPrecision
  use endfConstants
  use endfTable_class, only : endfTable, endf_bin_integral
  use pFUnit_mod

  implicit none

contains

  !!
  !! Test ENDF integration over a bin
  !!
@Test
  subroutine testEndfBinIntegral()
    real(defReal) :: x1, x0, y1, y0, val
    real(defReal), parameter :: TOL = 1.0E-5

    !! *** Test histogram
    x1 = 1.0_defReal
    x0 = 0.0_defReal
    y1 = 15.0_defReal
    y0 = 3.3_defReal
    ! @0.3
    val = 3.3_defReal * 0.3_defReal
    @assertEqual(val, endf_bin_integral(x0, x1, y0, y1, 0.3_defReal, histogramInterpolation), TOL*val)

    ! @1.0
    val = 3.3_defReal
    @assertEqual(val, endf_bin_integral(x0, x1, y0, y1, 1.0_defReal, histogramInterpolation), TOL*val)

    !! *** Test lin-lin Interpolation
    x1 = 2.0_defReal
    x0 = 1.0_defReal
    y1 = -1.0_defReal
    y0 = 1.0_defReal
    ! @1.5
    val = 0.25_defReal
    @assertEqual(val, endf_bin_integral(x0, x1, y0, y1, 1.5_defReal, linLinInterpolation), TOL*val)

    ! @2.0
    val = 0.0_defReal
    @assertEqual(val, endf_bin_integral(x0, x1, y0, y1, 2.0_defReal, linLinInterpolation), TOL)

    ! *** Test log-lin Interpolation (log in x-axis)
    ! y = 2*log(x)
    x1 = 5.0_defReal
    x0 = 1.0_defReal
    y1 = TWO * log(x1)
    y0 = 0.0_defReal
    ! @1.5
    val = 0.216395_defReal
    @assertEqual(val, endf_bin_integral(x0, x1, y0, y1, 1.5_defReal, logLinInterpolation), TOL*val)

    ! @5.0
    val = 8.09438_defReal
    @assertEqual(val, endf_bin_integral(x0, x1, y0, y1, 5.0_defReal, logLinInterpolation), TOL*val)

    ! *** Test lin-log Interpolation (lin in x-axis)
    ! y = 2*exp(x)
    x1 = 5.0_defReal
    x0 = 1.0_defReal
    y1 = TWO * exp(x1)
    y0 = TWO * exp(x0)
    ! @2.5
    val = TWO * (exp(2.5_defReal) - exp(x0))
    @assertEqual(val, endf_bin_integral(x0, x1, y0, y1, 2.5_defReal, linLogInterpolation), TOL*val)

    ! @5.0
    val = TWO * (exp(5.0_defReal) - exp(x0))
    @assertEqual(val, endf_bin_integral(x0, x1, y0, y1, 5.0_defReal, linLogInterpolation), TOL*val)


    ! *** Test log-log Interpolation
    ! Slope /= -1
    ! y = x^7
    x1 = 5.0_defReal
    x0 = 1.0_defReal
    y1 = x1**7
    y0 = x0**7
    ! @3.0
    val = (3.0_defReal**8 - x0**8) / 8.0_defReal
    @assertEqual(val, endf_bin_integral(x0, x1, y0, y1, 3.0_defReal, logLogInterpolation), TOL*val)

    ! @5.0
    val = (5.0_defReal**8 - x0**8) / 8.0_defReal
    @assertEqual(val, endf_bin_integral(x0, x1, y0, y1, 5.0_defReal, logLogInterpolation), TOL*val)

    ! Slope == -1
    ! y = 1/x
    x1 = 5.0_defReal
    x0 = 1.0_defReal
    y1 = ONE / x1
    y0 = ONE / x0
    ! @2.5
    val = log(2.5_defReal/x0)
    @assertEqual(val, endf_bin_integral(x0, x1, y0, y1, 2.5_defReal, logLogInterpolation), TOL*val)

    ! @5.0
    val = log(5.0_defReal/x0)
    @assertEqual(val, endf_bin_integral(x0, x1, y0, y1, 5.0_defReal, logLogInterpolation), TOL*val)
  end subroutine testEndfBinIntegral


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
