module tabularPdf_test

  use numPrecision
  use endfConstants
  use tabularPdf_class, only : tabularPdf
  use pFUnit_mod

  implicit none

@testCase
  type, extends(TestCase) :: test_tabularPdf
    private
    type(tabularPdf) :: pdf_lin
    type(tabularPdf) :: pdf_lin_CDF
    type(tabularPdf) :: pdf_hist
    type(tabularPdf) :: pdf_hist_CDF
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_tabularPdf


contains

  !!
  !! Sets up test_tabularPdf object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_tabularPdf), intent(inout) :: this
    integer(shortInt),parameter :: R = defReal
    real(defReal),dimension(4), parameter :: Grid = [ 1.0_R, 2.0_R, 3.0_R, 4.0_R]
    real(defReal),dimension(4), parameter :: PDF_L  = [ 0.0_R, 0.5_R, 0.5_R, 0.0_R]
    real(defReal),dimension(4), parameter :: CDF_L  = [ 0.0_R, 0.25_R, 0.75_R, 1.0_R]
    real(defReal),dimension(4), parameter :: PDF_H  = [ 0.3_R, 0.5_R, 0.2_R, 0.2_R]
    real(defReal),dimension(4), parameter :: CDF_H  = [ 0.0_R, 0.3_R, 0.8_R, 1.0_R]

    ! Initialise Linear interpolation table with and without provided CDF
    ! Symmetric trapezoidal PDF between 1-4
    call this % pdf_lin % init(Grid, PDF_L, tabPdfLinLin)
    call this % pdf_lin_CDF % init(Grid, PDF_L, CDF_L, tabPdfLinLin)

    ! Initialise histogram interpolation table with and without provided CDF
    ! Assymetric histogram with 3 bins (0.3, 0.5, 0.2) between 1-4
    call this % pdf_hist % init(Grid, PDF_H, tabPdfHistogram)
    call this % pdf_hist_CDF % init(Grid, PDF_H, CDF_H, tabPdfHistogram)

  end subroutine setUp

  !!
  !! Kills test_tabularPdf object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_tabularPdf), intent(inout) :: this

    call this % pdf_lin % kill()
    call this % pdf_lin_CDF % kill()
    call this % pdf_hist % kill()
    call this % pdf_hist_CDF % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  !!
  !!
  !!
@Test
  subroutine testReadingTable(this)
    class(test_tabularPdf), intent(inout) :: this
    real(defReal),parameter               :: TOL = 1.0E-9

    ! Test linear pdf
    @assertEqual(0.5_defReal, this % pdf_lin % probabilityOf(2.1_defReal), TOL)
    @assertEqual(0.5_defReal, this % pdf_lin_CDF % probabilityOf(2.1_defReal), TOL)

    @assertEqual(0.35_defReal, this % pdf_lin % probabilityOf(1.7_defReal), TOL)
    @assertEqual(0.35_defReal, this % pdf_lin_CDF % probabilityOf(1.7_defReal), TOL)

    @assertEqual(0.5_defReal, this % pdf_lin % probabilityOf(3.0_defReal), TOL)
    @assertEqual(0.5_defReal, this % pdf_lin_CDF % probabilityOf(3.0_defReal), TOL)

    ! Test Histogram PDF
    @assertEqual(0.5_defReal, this % pdf_hist % probabilityOf(2.1_defReal), TOL)
    @assertEqual(0.5_defReal, this % pdf_hist_CDF % probabilityOf(2.1_defReal), TOL)

    @assertEqual(0.3_defReal, this % pdf_hist % probabilityOf(1.7_defReal), TOL)
    @assertEqual(0.3_defReal, this % pdf_hist_CDF % probabilityOf(1.7_defReal), TOL)

    @assertEqual(0.2_defReal, this % pdf_hist % probabilityOf(3.01_defReal), TOL)
    @assertEqual(0.2_defReal, this % pdf_hist_CDF % probabilityOf(3.01_defReal), TOL)

  end subroutine testReadingTable

  !!
  !! Test getting bounds
  !!
@Test
  subroutine testGettingBounds(this)
    class(test_tabularPdf), intent(inout) :: this
    real(defReal),parameter               :: TOL = 1.0E-9
    real(defReal)                         :: bottom, top

    call this % pdf_lin % bounds(bottom, top)
    @assertEqual(1.0_defReal, bottom, TOL)
    @assertEqual(4.0_defReal, top, TOL)

    call this % pdf_lin_CDF % bounds(bottom, top)
    @assertEqual(1.0_defReal, bottom, TOL)
    @assertEqual(4.0_defReal, top, TOL)

    call this % pdf_hist % bounds(bottom, top)
    @assertEqual(1.0_defReal, bottom, TOL)
    @assertEqual(4.0_defReal, top, TOL)

    call this % pdf_hist_CDF % bounds(bottom, top)
    @assertEqual(1.0_defReal, bottom, TOL)
    @assertEqual(4.0_defReal, top, TOL)

  end subroutine testGettingBounds

  !!
  !! Test inversion of CDF
  !!
@Test
  subroutine testSample(this)
    class(test_tabularPdf), intent(inout) :: this
    real(defReal),parameter               :: TOL = 1.0E-9

    ! Linear PDFs
    @assertEqual(2.5_defReal, this % pdf_lin % sample(0.5_defReal), TOL)
    @assertEqual(2.5_defReal, this % pdf_lin_CDF % sample(0.5_defReal), TOL)

    @assertEqual(2.9_defReal, this % pdf_lin % sample(0.7_defReal), TOL)
    @assertEqual(2.9_defReal, this % pdf_lin_CDF % sample(0.7_defReal), TOL)

    @assertEqual(3.367544468_defReal, this % pdf_lin % sample(0.9_defReal), TOL)
    @assertEqual(3.367544468_defReal, this % pdf_lin_CDF % sample(0.9_defReal), TOL)

    @assertEqual(2.0_defReal, this % pdf_lin % sample(0.25_defReal), TOL)
    @assertEqual(2.0_defReal, this % pdf_lin_CDF % sample(0.25_defReal), TOL)

    @assertEqual(1.0_defReal, this % pdf_lin % sample(0.0_defReal), TOL)
    @assertEqual(1.0_defReal, this % pdf_lin_CDF % sample(0.0_defReal), TOL)

    @assertEqual(4.0_defReal, this % pdf_lin % sample(1.0_defReal), TOL)
    @assertEqual(4.0_defReal, this % pdf_lin_CDF % sample(1.0_defReal), TOL)

    ! Histograms PDFs
    @assertEqual(2.4_defReal, this % pdf_hist % sample(0.5_defReal), TOL)
    @assertEqual(2.4_defReal, this % pdf_hist_CDF % sample(0.5_defReal), TOL)

    @assertEqual(3.0_defReal, this % pdf_hist % sample(0.8_defReal), TOL)
    @assertEqual(3.0_defReal, this % pdf_hist_CDF % sample(0.8_defReal), TOL)

    @assertEqual(1.0_defReal, this % pdf_hist % sample(0.0_defReal), TOL)
    @assertEqual(1.0_defReal, this % pdf_hist_CDF % sample(0.0_defReal), TOL)

    @assertEqual(4.0_defReal, this % pdf_hist % sample(1.0_defReal), TOL)
    @assertEqual(4.0_defReal, this % pdf_hist_CDF % sample(1.0_defReal), TOL)

  end subroutine testSample

  !!
  !! Test inversion of CDF with bin
  !!
@Test
  subroutine testSampleWithBin(this)
    class(test_tabularPdf), intent(inout) :: this
    real(defReal),parameter               :: TOL = 1.0E-9
    integer(shortInt)                     :: bin

    ! Linear PDFs
    @assertEqual(2.5_defReal, this % pdf_lin % sample(0.5_defReal, bin), TOL)
    @assertEqual(2, bin)

    @assertEqual(3.367544468_defReal, this % pdf_lin_CDF % sample(0.9_defReal, bin), TOL)
    @assertEqual(3, bin)

    ! Histograms PDFs
    @assertEqual(2.4_defReal, this % pdf_hist % sample(0.5_defReal, bin), TOL)
    @assertEqual(2, bin)

    @assertEqual(4.0_defReal, this % pdf_hist_CDF % sample(1.0_defReal, bin), TOL)
    @assertEqual(3, bin)

  end subroutine testSampleWithBin


end module tabularPdf_test
