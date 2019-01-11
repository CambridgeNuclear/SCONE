module scoreMemory_test
  use numPrecision
  use scoreMemory_class, only : scoreMemory
  use pFUnit_mod

  implicit none

@testCase
  type, extends(TestCase) :: test_scoreMemory
    private
    type(scoreMemory) :: scoreMem
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_scoreMemory

contains

  !!
  !! Sets up test_scoreMemory object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_scoreMemory), intent(inout) :: this

    ! Allocate scoreMemory
    call this % scoreMem % init(5,1)

  end subroutine setUp

  !!
  !! Kills test_scoreMemory object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_scoreMemory), intent(inout) :: this

    call this % scoreMem % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Use test, Check if the mean and STD of the scores[1-10] is correct
  !!
@Test
  subroutine testScoring(this)
    class(test_scoreMemory), intent(inout) :: this
    real(defReal)                          :: mean,mean2, STD
    integer(shortInt)                      :: i,j

    do j=1,20   ! Loop over batches
      do i=1,10 ! Scoring Loop
        call this % scoreMem % score(real(j+i,defReal),1)
        call this % scoreMem % score(real(j+i,defReal),3)
      end do
      call this % scoreMem % normalise(ONE/12)
      call this % scoreMem % closeBatch(ONE)
    end do


    ! Test invalid index
    call this % scoreMem % getResult(mean, STD, -7)
    call this % scoreMem % getResult(mean2,-7)
    @assertEqual(ZERO, mean, 1.0E-9_defReal)
    @assertEqual(ZERO, mean, 1.0E-9_defReal)
    @assertEqual(ZERO, STD,  1.0E-9_defReal)

    ! Test empty bin
    call this % scoreMem % getResult(mean, STD, 2)
    call this % scoreMem % getResult(mean2, 2)
    @assertEqual(ZERO, mean, 1.0E-9_defReal)
    @assertEqual(ZERO, mean2, 1.0E-9_defReal)
    @assertEqual(ZERO, STD,  1.0E-9_defReal )

    ! Test bin 1
    call this % scoreMem % getResult(mean, STD, 1)
    call this % scoreMem % getResult(mean2, 1)
    @assertEqual(40.0_defReal/3, mean, 1.0E-9_defReal)
    @assertEqual(40.0_defReal/3, mean2, 1.0E-9_defReal)
    @assertEqual(1.102396379610246_defReal, STD,   1.0E-9_defReal)

    ! Test bin 3
    call this % scoreMem % getResult(mean, STD, 3)
    call this % scoreMem % getResult(mean, 3)
    @assertEqual(40.0_defReal/3, mean, 1.0E-9_defReal)
    @assertEqual(40.0_defReal/3, mean2, 1.0E-9_defReal)
    @assertEqual(1.102396379610246_defReal, STD,   1.0E-9_defReal)

  end subroutine testScoring

  !!
  !! Check correct behaviour of accumulate function
  !!
@Test
  subroutine testAccumulation(this)
    class(test_scoreMemory), intent(inout) :: this
    real(defReal)                          :: mean,mean2, STD
    integer(shortInt)                      :: i,j

    ! Accumulate number of results on bin 4
    call this % scoreMem % accumulate(1.1_defReal, 4)
    call this % scoreMem % accumulate(0.9_defReal, 4)
    call this % scoreMem % accumulate(1.15_defReal,4)

    ! Do some scoring
    do j=1,20   ! Loop over batches
      do i=1,10 ! Scoring Loop
        call this % scoreMem % score(real(j+i,defReal),1)
        call this % scoreMem % score(real(j+i,defReal),3)
      end do
      call this % scoreMem % normalise(ONE/12)
      call this % scoreMem % closeBatch(ONE)
    end do

    ! Do some more accumulation
    call this % scoreMem % accumulate(1.3_defReal, 4)
    call this % scoreMem % accumulate(0.95_defReal, 4)
    call this % scoreMem % accumulate(1.1_defReal,4)
    call this % scoreMem % accumulate(1_shortInt,4) ! Short Int
    call this % scoreMem % accumulate(1_longInt,4)  ! Long Int

    ! Get result
    call this % scoreMem % getResult(mean, STD, 4, 8)
    call this % scoreMem % getResult(mean2, 4, 8)

    ! Check accumulated results
    @assertEqual(1.0625_defReal, mean, 1.0E-9_defReal)
    @assertEqual(1.0625_defReal, mean2, 1.0E-9_defReal)
    @assertEqual(0.045069390943300_defReal, STD, 1.0E-9_defReal)

  end subroutine testAccumulation

  !!
  !! Test scoring with integers
  !!
@Test
  subroutine scoreIntegers(this)
    class(test_scoreMemory), intent(inout) :: this
    real(defReal)      :: mean, STD
    integer(shortInt)  :: i,j

    ! Do some scoring
    do j=1,20   ! Loop over batches
      do i=1,10 ! Scoring Loop
        call this % scoreMem % score(j+i,1)
        call this % scoreMem % score(int(j+i,longInt),3)
      end do
      call this % scoreMem % normalise(ONE/12)
      call this % scoreMem % closeBatch(ONE)
    end do

    ! Get result
    call this % scoreMem % getResult(mean, STD, 1)

    ! Check accumulated results
    ! Test bin 1
    call this % scoreMem % getResult(mean, STD, 1)
    @assertEqual(40.0_defReal/3, mean, 1.0E-9_defReal)
    @assertEqual(1.102396379610246_defReal, STD,   1.0E-9_defReal)

    ! Test bin 3
    call this % scoreMem % getResult(mean, STD, 3)
    @assertEqual(40.0_defReal/3, mean, 1.0E-9_defReal)
    @assertEqual(1.102396379610246_defReal, STD,   1.0E-9_defReal)

  end subroutine scoreIntegers
end module scoreMemory_test
