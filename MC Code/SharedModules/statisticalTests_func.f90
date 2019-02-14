module statisticalTests_func

  use numPrecision
  use genericProcedures, only : quickSort, fatalError

  implicit none


contains

  !!
  !! Computes Kolmogorov-Smirnov distance of CDFs for two samples of points
  !! 2nd sample is weighted by some weights
  !! Returns the estimate of p value
  !!
  !! Becouse some samples are weighthed an effective sample size
  !! is used for the weighted samples using Kish's formula taken from:
  !!   https://docs.displayr.com/wiki/Design_Effects_and_Effective_Sample_Size
  !!
  function twoSampleKS(sample1, sample2, wgts, D_ret) result(p)
    real(defReal),dimension(:),intent(inout) :: sample1
    real(defReal),dimension(:),intent(inout) :: sample2
    real(defReal),dimension(:),intent(inout) :: wgts
    real(defReal),intent(out),optional       :: D_ret
    real(defReal)                            :: p
    real(defReal)                            :: step1, step2, CDF1, CDF2, D
    real(defReal)                            :: N1, N2
    integer(shortInt)                        :: i1, i2, Top1, Top2
    logical(defBool)                         :: takeFrom2
    character(100), parameter :: Here = 'twoSampleKS (statisticalTests_func.f90)'

    ! Check size of sample2 and wgts
    if(size(sample2) /= size(wgts)) then
      call fatalError(Here,'Size of 2nd sample vector and vector of weights is diffrent')
    end if

    ! Catch input with an empty array of samples
    if(size(sample1) == 0 .or. size(sample2) == 0) then
      D = 1
      return
    end if

    ! Sort sample 1
    call quickSort(sample1)

    ! Sort sample 2 together with weights
    call quickSort(sample2, wgts)

    ! Calculate step size for samples array 1 and 2
    step1 = ONE / size(sample1)
    step2 = ONE / sum(wgts)

    ! Save size (top of samples1 and samples2)
    Top1 = size(sample1)
    Top2 = size(sample2)

    ! Initialise work variables
    CDF1 = ZERO
    CDF2 = ZERO
    i1 = 1
    i2 = 1
    D = ZERO

    ! Note that i1 & i2 are next bin indicaters, thus we loop until they reach Top + 1
    do while(i1 /= Top1 + 1 .and. i2 /= Top2 + 2)

      ! Calculate new maximum distance between CDFs
      D = max(D, abs(CDF1 - CDF2))

      ! Determine whether to take next value from samples2
      takeFrom2 = (sample2(i2) < sample1(i1) .and. i2 /= Top2 + 1) .or. i1 == Top1 + 1

      if(takeFrom2) then
        CDF2 = CDF2 + step2 * wgts(i2)
        i2 = i2 + 1

      else ! Take from sample1
        CDF1 = CDF1 + step1
        i1 = i1 + 1
      end if

    end do

    ! Calculate effective polulations accounting for weights
    N1 = real(Top1, defReal)
    N2 = sum(wgts) * sum(wgts) / sum(wgts * wgts)

    ! Calculate p-value
    p = exp(-TWO * D * D * N1 * N2/ (N1 + N2 ) )

    if(present(D_ret)) D_ret = D


  end function twoSampleKS



end module statisticalTests_func
