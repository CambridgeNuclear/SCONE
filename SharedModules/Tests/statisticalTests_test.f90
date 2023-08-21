module statisticalTests_test

  use numPrecision
  use statisticalTests_func, only : twoSampleKS
  use funit

  implicit none

contains

  !!
  !! Test Kolmogorov Smirnov Test on a simple case
  !! Taken from: http://influentialpoints.com/Training/kolmogorov-smirnov_test.htm
  !!
@Test
  subroutine testTwoSamplesKS()
    real(defReal), dimension(11) :: S1 = [51.0, 71.0, 42.0, 37.0, 51.0, 78.0, 51.0, 49.0, 56.0, 47.0, 58.0]
    real(defReal), dimension(5)  :: S2 = [45.0, 87.0, 123.0, 120.0, 70.0]
    real(defReal), dimension(5)  :: W  = ONE
    real(defReal)                :: p, D

    p = twoSampleKS(S1, S2, W, D)

    @assertEqual(0.6182_defReal, D, 0.00005_defReal)
    @assertEqual(0.07226_defReal, p, 0.00005_defReal )

  end subroutine testTwoSamplesKS

  !!
  !! Test Kolmogorov Smirnov Test on a simple case
  !! Modified for weights -> double weight of first element of S2
  !! Taken from: http://influentialpoints.com/Training/kolmogorov-smirnov_test.htm
  !!
@Test
  subroutine testTwoSamplesKS_withWgt()
    real(defReal), dimension(11) :: S1 = [51.0, 71.0, 42.0, 37.0, 51.0, 78.0, 51.0, 49.0, 56.0, 47.0, 58.0]
    real(defReal), dimension(5)  :: S2 = [45.0, 87.0, 123.0, 120.0, 70.0]
    real(defReal), dimension(5)  :: W  = ONE
    real(defReal)                :: p, D

    W(1) = TWO
    p = twoSampleKS(S1, S2, W, D)

    @assertEqual(0.48485_defReal, D, 0.00005_defReal)
    @assertEqual(0.22280_defReal, p, 0.00005_defReal )

  end subroutine testTwoSamplesKS_withWgt


end module statisticalTests_test
