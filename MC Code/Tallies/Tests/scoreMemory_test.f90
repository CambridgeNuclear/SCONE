module scoreMemory_test
  use numPrecision
  use genericProcedures, only : numToChar
  use scoreMemory_class, only : scoreMemory
  use pFUnit_mod

  implicit none

@testParameter(constructor = new_testNumber)
  type, extends(AbstractTestParameter) :: testNumber
    integer(shortInt) :: i
  contains
    procedure :: toString
  end type testNumber

@testCase(constructor=newTest)
  type, extends(ParameterizedTestCase) :: test_scoreMemory
    private
    integer(longInt)                            :: Ncycles
    integer(shortInt)                           :: batchSize
    real(defReal),dimension(:), allocatable     :: scores
    integer(shortInt), dimension(:),allocatable :: scoresInt

  end type test_scoreMemory


contains

  !!
  !! Build new test parameter form integer
  !!
  function new_testNumber(i) result (tstNum)
    integer(shortInt) :: i
    type(testNumber)  :: tstNum

    tstNum % i = i

  end function new_testNumber

  !!
  !! Write test parameter to string
  !!
  function toString(this) result(string)
    class(testNumber), intent(in) :: this
    character(:), allocatable :: string
    character(nameLen)        :: str

    write (str,*) this % i
    string = str

  end function toString

  !!
  !! Construct test case
  !!
  !!
  !!
  function newTest(testParam) result(tst)
    type(testNumber), intent(in)     :: testParam
    type(test_scoreMemory)           :: tst
    real(defReal),dimension(200)     :: random
    integer(shortInt)                :: seed, i
    integer(shortInt),parameter      :: A = 2469   ! Multiplier of LC PRNG
    integer(shortInt),parameter      :: M = 65521  ! Modulus of PRNG

    ! Load batchSize
    tst % batchSize = testParam % i
    tst % Ncycles   = 10 * tst % batchSize

    ! Generate a vector of 20 pseudo-random numbers in <0;1>
    ! Generator is not sophisticated but rebust
    seed = 9294
    do i=1,200
      seed = mod(A * seed , M)
      random(i)    = seed / real(M,defReal)
    end do

    ! Generate some scores and calculate their sum and sum of squares
    tst  % scores    = TWO + sin(PI * random - PI/2)
    tst % scoresInt = int(random * 100, shortInt)

  end function newTest

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test acoring for a case with batchSize == 1
  !! Look at the end of the file to find MATLAB SCRIPT used to generate reference values
  !!
@Test(cases=[1])
  subroutine testScoring(this)
    class(test_scoreMemory), intent(inout) :: this
    type(scoreMemory)                      :: mem
    integer(shortInt)                      :: i, j
    real(defReal)                          :: res1, res2, STD
    real(defReal), parameter :: TOL = 1.0E-9

    ! Initialise score memory
    call mem % init(7_longInt, 1, batchSize = this % batchSize)

    ! Test getting batchSize
    @assertEqual(this % batchSize, mem % getBatchSize(),'Test getBatchSize() :')

    ! Score in
    do i=1,10
      ! Score
      do j=20*(i-1)+1,20 * i
        call mem % score(this % scores(j), 1_longInt)
        call mem % score(this % scoresInt(j), 2_longInt)
        call mem % score(int(this % scoresInt(j),longInt),3_longInt)
        call mem % accumulate(this % scores(j), 4_longInt)
        call mem % accumulate(this % scoresInt(j), 5_longInt)
        call mem % accumulate(int(this % scoresInt(j),longInt),6_longInt)

      end do
      ! Close Cycle
      call mem % closeCycle(0.7_defReal)

    end do

    ! Get results from bin 1
    call mem % getResult(res1, 1_longInt)
    call mem % getResult(res2, STD, 1_longInt)

    @assertEqual(26.401471259728442_defReal, res1, TOL)
    @assertEqual(26.401471259728442_defReal, res2, TOL)
    @assertEqual(0.645969443981583_defReal, STD, TOL)

    ! Get results from bin 2
    call mem % getResult(res1, 2_longInt)
    call mem % getResult(res2, STD, 2_longInt)

    @assertEqual(623.0_defReal, res1, TOL)
    @assertEqual(623.0_defReal, res2, TOL)
    @assertEqual(27.982494527829360_defReal, STD, TOL)

    ! Get results from bin 3
    call mem % getResult(res1, 3_longInt)
    call mem % getResult(res2, STD, 3_longInt)

    @assertEqual(623.0_defReal, res1, TOL)
    @assertEqual(623.0_defReal, res2, TOL)
    @assertEqual(27.982494527829360_defReal, STD, TOL)

    ! Get results from bin 4
    call mem % getResult(res1, 4_longInt, 200)
    call mem % getResult(res2, STD, 4_longInt, 200)

    @assertEqual(1.885819375694888_defReal, res1, TOL)
    @assertEqual(1.885819375694888_defReal, res2, TOL)
    @assertEqual(0.049102082638055_defReal, STD, TOL)

    ! Get results from bin 5
    call mem % getResult(res1, 5_longInt, 200)
    call mem % getResult(res2, STD, 5_longInt, 200)

    @assertEqual(44.500000000000000_defReal, res1, TOL)
    @assertEqual(44.500000000000000_defReal, res2, TOL)
    @assertEqual(2.015580019267494_defReal, STD, TOL)

    ! Get results from bin 6
    call mem % getResult(res1, 6_longInt, 200)
    call mem % getResult(res2, STD, 6_longInt, 200)

    @assertEqual(44.500000000000000_defReal, res1, TOL)
    @assertEqual(44.500000000000000_defReal, res2, TOL)
    @assertEqual(2.015580019267494_defReal, STD, TOL)

    ! Get results from an empty bin 7
    call mem % getResult(res1, 7_longInt)
    call mem % getResult(res2, STD, 7_longInt)

    @assertEqual(ZERO, res1, TOL)
    @assertEqual(ZERO, res2, TOL)
    @assertEqual(ZERO, STD, TOL)

    ! Get results from invalid bins
    call mem % getResult(res1, -7_longInt)
    call mem % getResult(res2, STD, -7_longInt)

    @assertEqual(ZERO, res1, TOL)
    @assertEqual(ZERO, res2, TOL)
    @assertEqual(ZERO, STD, TOL)

    call mem % getResult(res1, 8_longInt)
    call mem % getResult(res2, STD, 8_longInt)

    @assertEqual(ZERO, res1, TOL)
    @assertEqual(ZERO, res2, TOL)
    @assertEqual(ZERO, STD, TOL)

    ! Free memor y
    call mem % kill()

  end subroutine testScoring

  !!
  !! Test lastCycle
  !! Ignors test parametrisation
  !!
@Test(cases=[1])
  subroutine testLastCycle(this)
    class(test_scoreMemory), intent(inout) :: this
    type(scoreMemory)                      :: mem
    integer(shortInt)                      :: i

    call mem % init(1_longInt, 1, batchSize = 8)

    ! Test getting batchSize
    @assertEqual(8, mem % getBatchSize(),'Test getBatchSize() :')

    do i=1,16
      if(i == 8 .or. i == 16) then
        @assertTrue( mem % lastCycle(), 'In cycle num: '//numToChar(i))
      else
        @assertFalse( mem % lastCycle(), 'In cycle num: '//numToChar(i))
      end if
      call mem % closeCycle(ONE)
    end do

    call mem % kill()

  end subroutine testLastCycle

  !!
  !! Test get score
  !! Ignore test parametrisation
  !!
@Test(cases=[1])
  subroutine testGetScore(this)
    class(test_scoreMemory), intent(inout) :: this
    type(scoreMemory)                      :: mem
    real(defReal),parameter :: TOL = 1.0E-9

    call mem % init(1_longInt, 1)

    call mem % score(ONE,1_longInt)
    call mem % score(ONE,1_longInt)
    call mem % score(ONE,1_longInt)

    @assertEqual(3*ONE, mem % getScore(1_longInt), TOL, 'Test getScore, valid bin:')
    @assertEqual(ZERO, mem % getScore(0_longInt), TOL, 'Test getScore, not +ve bin:')
    @assertEqual(ZERO, mem % getScore(2_longInt), TOL, 'Test getScore, too large bin:')

  end subroutine testGetScore

end module scoreMemory_test
!! MATLAB SCRIPT USED TO GENERATE REFERENCE VALUES
!clear
!rand = zeros(20,1);
!seed = 9294;
!
!%LCG Params
!A = 2469;
!M = 65521;
!
!for i=1:1:200
!  seed = mod(A * seed, M);
!  rand(i) = seed/M;
!end
!
!% Calculate scores vector
!scores = 2.0 + sin(pi() .* rand - pi()/2);
!scoresInt = floor(100.*rand);
!
!% Accumulate results
!resAcc = mean(scores)
!stdAcc = sqrt(var(scores)./200)
!
!resAccInt = mean(scoresInt)
!stdAccInt = sqrt(var(scoresInt)./200)
!
!% Reshape scores
!scores = reshape(scores,[20,10]);
!scores = sum(scores,1)* 0.7;
!res = mean(scores)
!std = sqrt(var(scores)./10)
!
!% Reshape scores
!scoresInt = reshape(scoresInt,[20,10]);
!scoresInt = sum(scoresInt,1)* 0.7;
!resInt = mean(scoresInt)
!stdInt = sqrt(var(scoresInt)./10)
