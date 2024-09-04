module RNG_test

  use numPrecision
  use RNG_class, only : RNG
  use pFUnit_mod

  implicit none

contains

  !!
  !! Test random numbers
  !!
@Test
  subroutine testRN()
    integer(shortInt), parameter :: N = 1000
    type(RNG)                    :: pRNG
    real(defReal), dimension(N)  :: rand
    integer(shortInt)            :: i

    call pRNG % init(int(z'5c3a84c9', longInt))

    do i=1,N
      rand(i) = pRNG % get()
    end do

    ! Check correcness
    @assertGreaterThanOrEqual(ONE, rand)
    @assertLessThanOrEqual(ZERO, rand)

  end subroutine testRN

  !!
  !! Test skip forward and backwards
  !!
@Test
  subroutine testSkip()
    type(RNG)         :: rand1
    type(RNG)         :: rand2
    integer(longInt)  :: seed
    real(defReal)     :: r_start, r2_start, r_end, r2_end
    integer(shortInt) :: i, N

    !! Initialise both RNGs to a nice number
    seed = int(z'5c3a84c9', longInt)
    call rand1 % init(seed)
    call rand2 % init(seed)

    !! Get initial random number
    r_start = rand1 % get()

    !! Move forward by 13456757 steps
    N = 13456757
    do i=1,N
      r_end = rand1 % get()
    end do

    ! Skip 2nd generator forward
    call rand2 % skip(int(N, longInt))
    r2_end = rand2 % get()

    ! Skip 2nd generator backwards. Must be 1 more becouse we drew a RN from generator
    call rand2 % skip(-int(N + 1, longInt))
    r2_start = rand2 % get()

    ! Verify values
    @assertEqual(r_end, r2_end)
    @assertEqual(r_start, r2_start)

  end subroutine testSkip



end module RNG_test
