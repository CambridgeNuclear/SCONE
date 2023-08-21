module sort_test
  use numPrecision
  use genericProcedures, only : swap, quickSort
  use funit

  implicit none


contains

  !!
  !! Test swaping
  !!
@Test
  subroutine testSwaping()
    real(defReal)     :: r1, r2, d1, d2
    integer(shortInt) :: i1, i2

    ! Swap diffrent integers
    i1 = -7
    i2 =  4

    call swap(i1, i2)

    @assertEqual(4_shortInt, i1)
    @assertEqual(-7_shortInt, i2)

    ! Swap same integers (If one tries XOR swap this can fail)
    i1 = 4
    i2 = 4
    call swap(i1, i2)
    @assertEqual(4_shortInt, i1)
    @assertEqual(4_shortInt, i2)

    ! Swap reals
    r1 = -7.3_defReal
    r2 = 2.3_defReal

    call swap(r1, r2)

    @assertEqual(2.3_defReal, r1)
    @assertEqual(-7.3_defReal, r2)

    ! Swap same reals
    r1 = 2.3_defReal
    r2 = 2.3_defReal

    call swap(r1, r2)

    @assertEqual(2.3_defReal, r1)
    @assertEqual(2.3_defReal, r2)

    ! Swap pair of reals
    r1 = 2.3_defReal
    r2 = 0.3_defReal
    d1 = 1.7_defReal
    d2 = 0.4_defReal

    call swap(r1, r2, d1, d2)

    @assertEqual(1.7_defReal, r1)
    @assertEqual(0.4_defReal, r2)
    @assertEqual(2.3_defReal, d1)
    @assertEqual(0.3_defReal, d2)

    ! Swap same pair of reals
    r1 = 2.3_defReal
    r2 = 0.3_defReal
    d1 = 2.3_defReal
    d2 = 0.3_defReal

    call swap(r1, r2, d1, d2)

    @assertEqual(2.3_defReal, r1)
    @assertEqual(0.3_defReal, r2)
    @assertEqual(2.3_defReal, d1)
    @assertEqual(0.3_defReal, d2)

  end subroutine testSwaping

  !!
  !! Test sorting integer array
  !!
@Test
  subroutine testQuickSortInt()
    integer(shortInt), dimension(4) :: I = [2, 1, 4, 3]
    integer(shortInt), dimension(:), allocatable :: IA

    ! Ordinary sort
    call quickSort(I)
    @assertEqual([1, 2, 3, 4], I)

    ! Sorted array
    I = [1, 2, 3, 4]
    call quickSort(I)
    @assertEqual([1, 2, 3, 4], I)

    ! Uniform Array
    I = [1, 1, 1, 1]
    call quickSort(I)
    @assertEqual([1, 1, 1, 1], I)

    ! 1 Element array
    IA = [1]
    call quickSort(IA)
    @assertEqual([1],IA)

    ! Empty array
    deallocate(IA)
    allocate(IA(0))

    call quickSort(IA)
    @assertEqual(0, size(IA))

  end subroutine testQuickSortInt

  !!
  !! Test sorting real array
  !!
@Test
  subroutine testQuickSortReal()
    real(defReal), dimension(4) :: R = [2.1_defReal, 1.4_defReal, 4.9_defReal, 3.4_defReal]
    real(defReal), dimension(:), allocatable :: RA

    ! Ordinary sort
    call quickSort(R)
    @assertEqual([1.4_defReal, 2.1_defReal, 3.4_defReal, 4.9_defReal], R)

    ! Sorted array
    R = [1.4_defReal, 2.1_defReal, 3.4_defReal, 4.9_defReal]
    call quickSort(R)
    @assertEqual([1.4_defReal, 2.1_defReal, 3.4_defReal, 4.9_defReal], R)

    ! Uniform Array
    R = [1.0_defReal, 1.0_defReal, 1.0_defReal, 1.0_defReal]
    call quickSort(R)
    @assertEqual([1.0_defReal, 1.0_defReal, 1.0_defReal, 1.0_defReal], R)

    ! 1 Element array
    RA = [1.0_defReal]
    call quickSort(RA)
    @assertEqual([1.0_defReal],RA)

    ! Empty array
    deallocate(RA)
    allocate(RA(0))

    call quickSort(RA)
    @assertEqual(0, size(RA))

  end subroutine testQuickSortReal

  !!
  !! Test sorting two arrays of reals by keys in first array
  !!
@Test
  subroutine testQuickSortRealReal()
    real(defReal), dimension(4) :: R1 = [2.1_defReal, 1.4_defReal, 4.9_defReal, 3.4_defReal]
    real(defReal), dimension(4) :: R2 = [1.0_defReal, 2.0_defReal, 3.0_defReal, 4.0_defReal]

    ! Ordinary sort
    call quickSort(R1, R2)
    @assertEqual([1.4_defReal, 2.1_defReal, 3.4_defReal, 4.9_defReal], R1)
    @assertEqual([2.0_defReal, 1.0_defReal, 4.0_defReal, 3.0_defReal], R2)

    ! Sorted array
    R1 = [1.4_defReal, 2.1_defReal, 3.4_defReal, 4.9_defReal]
    R2 = [1.0_defReal, 2.0_defReal, 3.0_defReal, 4.0_defReal]

    call quickSort(R1, R2)
    @assertEqual([1.4_defReal, 2.1_defReal, 3.4_defReal, 4.9_defReal], R1)
    @assertEqual([1.0_defReal, 2.0_defReal, 3.0_defReal, 4.0_defReal], R2)

    ! Uniform array
    R1 = [1.0_defReal, 1.0_defReal, 1.0_defReal, 1.0_defReal]
    R2 = [1.0_defReal, 2.0_defReal, 3.0_defReal, 4.0_defReal]

    call quickSort(R1, R2)
    @assertEqual([1.0_defReal, 1.0_defReal, 1.0_defReal, 1.0_defReal], R1)
    @assertEqual([1.0_defReal, 2.0_defReal, 3.0_defReal, 4.0_defReal], R2)

  end subroutine testQuickSortRealReal


end module sort_test
