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

!! Parameters
  real(defReal), dimension(*), parameter :: scores1  = [-0.806853587_defReal,&
   -0.264065733_defReal, 0.748995766_defReal, 0.649558345_defReal, -1.729071673_defReal,&
   0.448515466_defReal, -0.557721502_defReal, 0.324925031_defReal, 1.379322451_defReal,&
   -0.247145475_defReal, -0.078875407_defReal, 1.078486121_defReal]

  real(defReal), dimension(*), parameter :: scores2  = [-0.034871774_defReal,&
   1.67602827_defReal, -0.128109861_defReal, 0.143275657_defReal, 0.99694274_defReal,&
   -1.411425667_defReal, -1.043824927_defReal, 0.438639436_defReal, -0.716936163_defReal,&
   0.888037098_defReal, -1.946112113_defReal, -0.762986825_defReal]

  real(defReal), dimension(*), parameter :: scores3  = [-0.217173463_defReal,&
  -1.100374936_defReal, 0.076404552_defReal, -0.41166647_defReal, 1.820644781_defReal,&
  0.117637644_defReal, -0.468246467_defReal, -1.049153503_defReal, -1.603418291_defReal,&
  1.083670697_defReal, -1.582420998_defReal, -0.581396211_defReal]

  real(defReal), dimension(*), parameter :: scores4  = [1.827634112_defReal,&
   1.860832525_defReal, 0.32182295_defReal, -1.005135626_defReal, -1.169059025_defReal,&
   1.626930992_defReal, -0.395332559_defReal, 0.437490287_defReal, -0.572598161_defReal,&
   0.939834964_defReal, -0.873795889_defReal, -1.487200675_defReal]

  real(defReal), dimension(*), parameter :: scores5  = [0.474213109_defReal,&
  1.239262389_defReal, -1.958608629_defReal, 1.747440442_defReal, -1.885318032_defReal,&
  -1.893646658_defReal, -1.940449615_defReal, 0.097153478_defReal, -0.137848924_defReal,&
  -0.44852551_defReal, -0.222650372_defReal, -1.996516242_defReal]

  real(defReal), dimension(*), parameter :: scores6  = [-1.585254618_defReal,&
  -0.103171792_defReal, -0.021466453_defReal, 0.880693845_defReal, 1.848642133_defReal,&
  -0.458288777_defReal, -1.543692062_defReal, -1.672100383_defReal, -0.150767892_defReal,&
  -1.882317987_defReal, -0.104484479_defReal, 1.333161081_defReal]

  real(defReal), dimension(*), parameter :: scores7  = [-0.494495838_defReal,&
  0.092782855_defReal, -1.939612951_defReal, 0.628355122_defReal, -0.868568289_defReal,&
  1.213031586_defReal, -0.479822236_defReal, -1.374732604_defReal, 0.869319694_defReal,&
  -1.571185737_defReal, 0.084051731_defReal, 1.317419474_defReal]

  real(defReal), dimension(*), parameter :: scores8  = [-1.297821846_defReal,&
  -0.490394877_defReal, 1.982176508_defReal, -0.821228721_defReal, -1.958716843_defReal,&
  -0.658591026_defReal, -0.162838986_defReal, 0.67008506_defReal, 1.99566651_defReal,&
  -0.443620463_defReal, -1.220572251_defReal, -0.387782805_defReal]

  real(defReal), dimension(*), parameter :: scores9  = [0.194196711_defReal,&
  1.533681665_defReal, 1.632933828_defReal, -1.717756433_defReal, 1.463262014_defReal,&
  -0.850636142_defReal, 0.409183378_defReal, 1.889650079_defReal, 0.30002653_defReal,&
  -0.190091317_defReal, 1.952906198_defReal, -0.602948161_defReal]

  real(defReal), dimension(*), parameter :: scores10 = [0.935316855_defReal,&
  -1.555582645_defReal, 0.066962677_defReal, 0.933512898_defReal, -1.606313364_defReal,&
  -1.216797735_defReal, 0.017590584_defReal, 1.836462286_defReal, 1.452971352_defReal,&
  0.459380287_defReal, -1.339129295_defReal, -1.340975458_defReal]

  real(defReal), parameter :: Sample_mean = -0.12002646_defReal
  real(defReal), parameter :: Sample_STD  = 0.093580078_defReal


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
    real(defReal)                          :: mean, STD
    integer(shortInt)                      :: i

    ! Score all batches on idx 1 and 3
    ! Batch 1
    do i=1,size(scores1)
      call this % scoreMem % score(scores1(i),1)
      call this % scoreMem % score(scores1(i),3)
    end do
    call this % scoreMem % normalise(ONE/12)
    call this % scoreMem % closeBatch(ONE)
    ! Batch 2
    do i=1,size(scores2)
      call this % scoreMem % score(scores2(i),1)
      call this % scoreMem % score(scores2(i),3)
    end do
    call this % scoreMem % normalise(ONE/12)
    call this % scoreMem % closeBatch(ONE)
    ! Batch 3
    do i=1,size(scores1)
      call this % scoreMem % score(scores3(i),1)
      call this % scoreMem % score(scores3(i),3)
    end do
    call this % scoreMem % normalise(ONE/12)
    call this % scoreMem % closeBatch(ONE)
    ! Batch 4
    do i=1,size(scores4)
      call this % scoreMem % score(scores4(i),1)
      call this % scoreMem % score(scores4(i),3)
    end do
    call this % scoreMem % normalise(ONE/12)
    call this % scoreMem % closeBatch(ONE)
    ! Batch 5
    do i=1,size(scores5)
      call this % scoreMem % score(scores5(i),1)
      call this % scoreMem % score(scores5(i),3)
    end do
    call this % scoreMem % normalise(ONE/12)
    call this % scoreMem % closeBatch(ONE)
    ! Batch 6
    do i=1,size(scores6)
      call this % scoreMem % score(scores6(i),1)
      call this % scoreMem % score(scores6(i),3)
    end do
    call this % scoreMem % normalise(ONE/12)
    call this % scoreMem % closeBatch(ONE)
    ! Batch 7
    do i=1,size(scores7)
      call this % scoreMem % score(scores7(i),1)
      call this % scoreMem % score(scores7(i),3)
    end do
    call this % scoreMem % normalise(ONE/12)
    call this % scoreMem % closeBatch(ONE)
    ! Batch 8
    do i=1,size(scores8)
      call this % scoreMem % score(scores8(i),1)
      call this % scoreMem % score(scores8(i),3)
    end do
    call this % scoreMem % normalise(ONE/12)
    call this % scoreMem % closeBatch(ONE)
    ! Batch 9
    do i=1,size(scores9)
      call this % scoreMem % score(scores9(i),1)
      call this % scoreMem % score(scores9(i),3)
    end do
    call this % scoreMem % normalise(ONE/12)
    call this % scoreMem % closeBatch(ONE)
    ! Batch 10
    do i=1,size(scores10)
      call this % scoreMem % score(scores10(i),1)
      call this % scoreMem % score(scores10(i),3)
    end do
    call this % scoreMem % normalise(ONE/12)
    call this % scoreMem % closeBatch(ONE)

    ! Test invalid index
    call this % scoreMem % getResult(mean, STD, -7)
    @assertEqual(ZERO, mean, 1.0E-9_defReal)
    @assertEqual(ZERO, STD,  1.0E-9_defReal)

    ! Test empty bin
    call this % scoreMem % getResult(mean, STD, 2)
    @assertEqual(ZERO, mean, 1.0E-9_defReal)
    @assertEqual(ZERO, STD,  1.0E-9_defReal )

    ! Test bin 1
    call this % scoreMem % getResult(mean, STD, 1)
    @assertEqual(Sample_mean, mean, 1.0E-9_defReal)
    @assertEqual(Sample_STD, STD,   1.0E-9_defReal)

    ! Test bin 3
    call this % scoreMem % getResult(mean, STD, 3)
    @assertEqual(Sample_mean, mean, 1.0E-9_defReal)
    @assertEqual(Sample_STD, STD,   1.0E-9_defReal)

  end subroutine testScoring

  !!
  !! Check correct behaviour of accumulate function
  !!
@Test
  subroutine testAccumulation(this)
    class(test_scoreMemory), intent(inout) :: this
    real(defReal)                          :: mean, STD
    integer(shortInt)                      :: i

    ! Accumulate number of results on bin 4
    call this % scoreMem % accumulate(1.1_defReal, 4)
    call this % scoreMem % accumulate(0.9_defReal, 4)
    call this % scoreMem % accumulate(1.15_defReal,4)

    ! Do some scoring
    ! Batch 7
    do i=1,size(scores7)
      call this % scoreMem % score(scores7(i),1)
      call this % scoreMem % score(scores7(i),3)
    end do
    call this % scoreMem % normalise(ONE/12)
    call this % scoreMem % closeBatch(ONE)
    ! Batch 8
    do i=1,size(scores8)
      call this % scoreMem % score(scores8(i),1)
      call this % scoreMem % score(scores8(i),3)
    end do
    call this % scoreMem % normalise(ONE/12)
    call this % scoreMem % closeBatch(ONE)
    ! Batch 9
    do i=1,size(scores9)
      call this % scoreMem % score(scores9(i),1)
      call this % scoreMem % score(scores9(i),3)
    end do
    call this % scoreMem % normalise(ONE/12)
    call this % scoreMem % closeBatch(ONE)

    ! Do some more accumulation
    call this % scoreMem % accumulate(1.3_defReal, 4)
    call this % scoreMem % accumulate(0.95_defReal, 4)
    call this % scoreMem % accumulate(1.1_defReal,4)
    call this % scoreMem % accumulate(1_shortInt,4) ! Short Int
    call this % scoreMem % accumulate(1_longInt,4)  ! Long Int

    ! Get result
    call this % scoreMem % getResult(mean, STD, 4, 8)

    ! Check accumulated results
    @assertEqual(1.0625_defReal, mean, 1.0E-9_defReal)
    @assertEqual(0.045069390943300_defReal, STD, 1.0E-9_defReal)

  end subroutine testAccumulation

  !!
  !! Test scoring with integers
  !!
@Test
  subroutine scoreIntegers(this)
    class(test_scoreMemory), intent(inout) :: this
    real(defReal) :: mean, STD

    ! Do some scoring
    call this % scoreMem % score(3_shortInt,1)
    call this % scoremem % score(4_longInt,1)
    call this % scoreMem % closeBatch(ONE)

    call this % scoreMem % score(7_shortInt,1)
    call this % scoremem % score(3_longInt,1)
    call this % scoreMem % closeBatch(ONE)

    call this % scoreMem % score(2_shortInt,1)
    call this % scoremem % score(6_longInt,1)
    call this % scoreMem % closeBatch(ONE)

    ! Get result
    call this % scoreMem % getResult(mean, STD, 1)

    ! Check accumulated results
    @assertEqual(8.333333333333334_defReal, mean, 1.0E-9_defReal)
    @assertEqual(0.881917103688197_defReal, STD, 1.0E-9_defReal)

  end subroutine scoreIntegers
end module scoreMemory_test
