module heapQueue_test
  use numPrecision
  use heapQueue_class, only: heapQueue
  use funit

  implicit none

contains

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test with simple sequence without reaching maximum size
  !!
@Test
  subroutine testBelowMaximum()
    type(heapQueue)   :: hq
    integer(shortInt) :: i
    real(defReal), dimension(*), parameter :: seq = [2.0_defReal, 3.0_defReal, 1.0_defReal, 4.0_defReal, 5.0_defReal]

    call hq % init(8)

    @assertEqual(hq % getSize(), 0)

    do i = 1, size(seq)
      call hq % pushReplace(seq(i))
    end do

    ! Check that the maximum value is the maximum value in the sequence
    @assertEqual(hq % maxValue(), maxval(seq))
    @assertEqual(hq % getSize(), size(seq))

  end subroutine testBelowMaximum

  !!
  !! Test the intended use case
  !!
  !! It is not a very good test but I had no better idea at the moment [MAK]
  !!
  @Test
  subroutine testAboveMaximum()
    type(heapQueue)   :: hq
    integer(shortInt) :: i
    real(defReal)     :: val
    real(defReal), dimension(*), parameter :: seq = [2.0_defReal, 3.0_defReal, 1.0_defReal, 1.4_defReal, 5.0_defReal]

    call hq % init(3)

    ! Push upper bound
    call hq % pushReplace(1000.0_defReal)

    do i = 1, size(seq)
      val = seq(i)
      if (val < hq % maxValue()) call hq % pushReplace(seq(i))
    end do

    ! Check that the threshold is correct
    @assertEqual(hq % maxValue(), 2.0_defReal)

  end subroutine testAboveMaximum

end module heapQueue_test
