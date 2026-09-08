module kalbachTable_test

  use endfConstants
  use funit
  use kalbachTable_class, only : kalbachTable
  use numPrecision

  implicit none

contains

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Loaded CDF ends inside the acceptance tolerance but below 1.0. Random 
  !! numbers in (cdf(N), 1.0] must still produce a valid sample instead of 
  !! a mid-run search failure. Covers the generator top edge rand = 1.0 as 
  !! well.
  !!
@Test
  subroutine testSampleAboveCDFEnd()
    real(defReal)                          :: A, R, x
    type(kalbachTable)                     :: tab
    real(defReal), dimension(3), parameter :: X_GRID = [1.0_defReal, 2.0_defReal, 3.0_defReal]
    real(defReal), dimension(3), parameter :: PDF    = [0.5_defReal, 0.5_defReal, 0.5_defReal]
    real(defReal), dimension(3), parameter :: CDF    = [0.0_defReal, 0.5_defReal, &
                                                        1.0_defReal - 1.0E-7_defReal]
    real(defReal), dimension(3), parameter :: R_GRID = [0.1_defReal, 0.2_defReal, 0.3_defReal]
    real(defReal), dimension(3), parameter :: A_GRID = [0.4_defReal, 0.5_defReal, 0.6_defReal]

    ! Initialise CDF.
    call tab % init(X_GRID, PDF, CDF, R_GRID, A_GRID, tabPdfLinLin)

    ! Sample using a random draw just below 1.0 and check that the number 
    ! sampled is valid.
    call tab % sample(1.0_defReal - 1.0E-8_defReal, x, R, A)
    @assertTrue(x == x)
    @assertTrue(x >= 1.0_defReal .and. x <= 3.0_defReal)

    ! Sample using a random draw of exactly 1.0 and check that the number 
    ! sampled is valid.
    call tab % sample(1.0_defReal, x, R, A)
    @assertTrue(x == x)
    @assertTrue(x >= 1.0_defReal .and. x <= 3.0_defReal)

    ! Clean up.
    call tab % kill()

  end subroutine testSampleAboveCDFEnd

  !!
  !! Loaded CDF begins inside the acceptance tolerance but above 0.0. Random
  !! numbers in [0.0, cdf(1)) must still produce a valid sample instead of a
  !! mid-run search failure.
  !!
@Test
  subroutine testSampleBelowCDFStart()
    real(defReal)                          :: A, R, x
    type(kalbachTable)                     :: tab
    real(defReal), dimension(3), parameter :: X_GRID = [1.0_defReal, 2.0_defReal, 3.0_defReal]
    real(defReal), dimension(3), parameter :: PDF    = [0.5_defReal, 0.5_defReal, 0.5_defReal]
    real(defReal), dimension(3), parameter :: CDF    = [1.0E-7_defReal, 0.5_defReal, 1.0_defReal]
    real(defReal), dimension(3), parameter :: R_GRID = [0.1_defReal, 0.2_defReal, 0.3_defReal]
    real(defReal), dimension(3), parameter :: A_GRID = [0.4_defReal, 0.5_defReal, 0.6_defReal]

    ! Initialise CDF.
    call tab % init(X_GRID, PDF, CDF, R_GRID, A_GRID, tabPdfLinLin)

    ! Sample using a random draw just above 0.0 and check that the number 
    ! sampled is valid.
    call tab % sample(1.0E-8_defReal, x, R, A)
    @assertTrue(x == x)
    @assertTrue(x >= 1.0_defReal .and. x <= 3.0_defReal)

    ! Clean up.
    call tab % kill()

  end subroutine testSampleBelowCDFStart

  !!
  !! linLin bin ending at zero density: the exact-arithmetic discriminant at the
  !! bin top is 0.0, so floating-point evaluation can round it below 0.0 and feed
  !! sqrt() a negative argument. The grid below makes the rounding land negative
  !! in IEEE double precision, so the sample must come from the fallback branch.
  !!
@Test
  subroutine testNegativeDiscriminant()
    real(defReal)                          :: A, R, x
    type(kalbachTable)                     :: tab
    real(defReal), dimension(2), parameter :: X_GRID = [1.0_defReal, 1.003_defReal]
    real(defReal), dimension(2), parameter :: PDF    = [2.0_defReal/0.003_defReal, 0.0_defReal]
    real(defReal), dimension(2), parameter :: CDF    = [0.0_defReal, 1.0_defReal]
    real(defReal), dimension(2), parameter :: R_GRID = [0.1_defReal, 0.2_defReal]
    real(defReal), dimension(2), parameter :: A_GRID = [0.4_defReal, 0.5_defReal]

    ! Initialise CDF.
    call tab % init(X_GRID, PDF, CDF, R_GRID, A_GRID, tabPdfLinLin)

    ! Sample using a random draw giving a negative discriminant and check that 
    ! the numbers sampled are valid.
    call tab % sample(1.0_defReal, x, R, A)
    @assertTrue(x == x)
    @assertTrue(x >= 1.0_defReal .and. x <= 1.003_defReal)
    @assertTrue(R == R)
    @assertTrue(A == A)

    ! Clean up.
    call tab % kill()

  end subroutine testNegativeDiscriminant

  !!
  !! Histogram distribution with a trailing zero-density bin: the computed CDF is
  !! flat over the last bin, so rand = 1.0 resolves onto the plateau where the
  !! unguarded inversion evaluates 0.0/0.0. The sample must return the plateau edge.
  !!
@Test
  subroutine testZeroDensityPlateau()
    real(defReal)                          :: A, R, x
    type(kalbachTable)                     :: tab
    real(defReal), parameter               :: TOL    = 1.0E-9_defReal
    real(defReal), dimension(3), parameter :: X_GRID = [1.0_defReal, 2.0_defReal, 3.0_defReal]
    real(defReal), dimension(3), parameter :: PDF    = [1.0_defReal, 0.0_defReal, 0.0_defReal]
    real(defReal), dimension(3), parameter :: R_GRID = [0.1_defReal, 0.2_defReal, 0.3_defReal]
    real(defReal), dimension(3), parameter :: A_GRID = [0.4_defReal, 0.5_defReal, 0.6_defReal]

    ! Initialise CDF.
    call tab % init(X_GRID, PDF, R_GRID, A_GRID, tabPdfHistogram)

    ! Sample using a random draw on the plateau and check that the numbers sampled
    ! correspond to those at the plateau edge.
    call tab % sample(1.0_defReal, x, R, A)
    @assertEqual(2.0_defReal, x, TOL)
    @assertEqual(0.2_defReal, R, TOL)
    @assertEqual(0.5_defReal, A, TOL)

    ! Clean up.
    call tab % kill()

  end subroutine testZeroDensityPlateau

end module kalbachTable_test