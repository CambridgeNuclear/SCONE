module endfTable_test
  use numPrecision
  use endfConstants
  use endfTable_class, only : endfTable, endf_bin_integral
  use funit

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

    ! Test reload
    call tab % reloadY([ONE, ONE, ONE, ONE])
    @assertEqual(ONE, tab % at(0.0_defReal), TOL)
    @assertEqual(ONE, tab % at(3.0_defReal), TOL)

  end subroutine testLinLinTable

  !!
  !! Test integration in lin-lin table
  !!
@Test
  subroutine testLinLinIntegral()
    real(defReal), dimension(*), parameter :: x = [1.0_defReal, 2.0_defReal, 2.5_defReal, 3.0_defReal]
    real(defReal), dimension(*), parameter :: y = [0.0_defReal, 1.0_defReal, 0.5_defReal, 0.0_defReal]
    type(endfTable)                        :: tab
    real(defReal), dimension(4)            :: out, val, points
    real(defReal),parameter :: TOL = 1.0E-6

    call tab % init(x, y)

    ! Integral at grid points
    out = tab % integral(x)
    val = [0.0_defReal, 0.5_defReal, 0.875_defReal, 1.0_defReal]
    @assertEqual(val, out, TOL)

    ! Integral of grid points
    points = [1.25_defReal, 1.75_defReal, 2.6_defReal, 2.9_defReal]
    val    = [3.125E-2_defReal, 0.28125_defReal, 0.92_defReal, 0.995_defReal]
    out = tab % integral(points)
    @assertEqual(val, out, TOL)


  end subroutine testLinLinIntegral


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
  !! Test integration with single (log-log) intrepolation egion
  !!
  !! Polynomial functions
  !!  x^2 if x < 1
  !!  x^(-1) if x >= 1
  !!
@Test
  subroutine testOneRegionIntegral()
    real(defReal), dimension(6), parameter :: x = [0.01_defReal, 0.5_defReal, 1.0_defReal,&
                                                   2.0_defReal, 4.0_defReal, 8.0_defReal]
    real(defReal), dimension(6), parameter :: y = [1.0E-4_defReal, 0.25_defReal, 1.0_defReal, &
                                                   0.5_defReal, 0.25_defReal, 0.125_defReal ]
    type(endfTable)                        :: tab
    real(defReal), dimension(6)            :: out, val, p
    real(defReal),parameter                :: TOL = 1.0E-6

    call tab % init(x, y, [6], [loglogInterpolation])

    ! Integral at grid points
    val(1) = ZERO
    val(2) = (x(2)**3 - x(1)**3) / 3.0_defReal
    val(3) = (ONE - x(1)**3) / 3.0_defReal
    val(4) = val(3) + log(x(4)/x(3))
    val(5) = val(3) + log(x(5)/x(3))
    val(6) = val(3) + log(x(6)/x(3))
    out = tab % integral(x)
    @assertEqual(val, out, TOL)

    ! Integral in between
    p = [0.1_defReal, 0.75_defReal, 0.9_defReal, 1.1_defReal, 1.4_defReal, 6.0_defReal]
    val(1) = (p(1)**3 - x(1)**3) / 3.0_defReal
    val(2) = (p(2)**3 - x(1)**3) / 3.0_defReal
    val(3) = (p(3)**3 - x(1)**3) / 3.0_defReal
    val(4) = ONE/3.0_defReal + log(p(4)/x(3)) - x(1)**3 / 3.0_defReal
    val(5) = ONE/3.0_defReal + log(p(5)/x(3)) - x(1)**3 / 3.0_defReal
    val(6) = ONE/3.0_defReal + log(p(6)/x(3)) - x(1)**3 / 3.0_defReal
    out = tab % integral(p)
    @assertEqual(val, out, TOL)

  end subroutine testOneRegionIntegral

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

  !!
  !! Test Integration with multiple interpolation regions
  !!  3 regions hist, lin and
  !!   Histogram for x < 2.0
  !!   Lin-Lin   for x < 3.0  y = x - 2
  !!   Log-Log   for x > 3.0  y = 3/x
  !!
  !!
@Test
  subroutine testMultiRegionIntegral()
    real(defReal), dimension(7), parameter :: x = [1.0_defReal, 2.0_defReal, 2.0_defReal, &
                                                   2.5_defReal, 3.0_defReal, 4.0_defReal, 7.0_defReal]
    real(defReal), dimension(7)            :: y
    type(endfTable)                        :: tab
    real(defReal), dimension(7)            :: out, val
    real(defReal), dimension(3)            :: p, val2
    real(defReal),parameter                :: TOL = 1.0E-6

    ! Calculate y values
    y(1) = 2.0_defReal
    y(2) = -4.0_defReal  ! Irrelevant for hist-interp
    y(3:5) = x(3:5) - 2.0_defReal
    y(6:7) = 3.0_defReal / x(6:7)

    ! Create table
    call tab % init(x, y, [2, 5, 7], [histogramInterpolation, linLinInterpolation, loglogInterpolation])

    ! Integral at grid points
    val(1) = ZERO
    val(2) = TWO
    val(3) = TWO
    val(4) = TWO + 0.125_defReal
    val(5) = TWO + HALF
    val(6) = TWO + HALF + 3*log(x(6)/x(5))
    val(7) = TWO + HALF + 3*log(x(7)/x(5))
    out = tab % integral(x)
    @assertEqual(val, out, TOL)

    ! Integral out grid points
    p = [1.5_defReal, 2.9_defReal, 3.2_defReal]
    val2(1) = ONE
    val2(2) = TWO + 0.405_defReal
    val2(3) = TWO + HALF + 3*log(p(3)/x(5))

  end subroutine testMultiRegionIntegral


end module endfTable_test
