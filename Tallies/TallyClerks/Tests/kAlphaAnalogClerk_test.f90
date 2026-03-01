module kAlphaAnalogClerk_test

  use numPrecision
  use tallyResult_class,       only : tallyResult
  use kAlphaAnalogClerk_class, only : kAlphaAnalogClerk, kAlphaResult
  use particle_class,          only : particle
  use particleDungeon_class,   only : particleDungeon
  use dictionary_class,        only : dictionary
  use scoreMemory_class,       only : scoreMemory
  use outputFile_class,        only : outputFile
  use funit

  implicit none

@testCase
  type, extends(TestCase) :: test_kAlphaAnalogClerk
    private
    type(kAlphaAnalogClerk) :: clerk
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_kAlphaAnalogClerk

contains

  !!
  !! Sets up test_kAlphaAnalogClerk object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_kAlphaAnalogClerk), intent(inout) :: this
    type(dictionary)                             :: dict
    character(nameLen)                           :: name

    call dict % init(2)
    call dict % store('alpha_0', ONE)
    call this % clerk % init(dict, name)

    call dict % kill()

  end subroutine setUp

  !!
  !! Kills test_kAlphaAnalogClerk object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_kAlphaAnalogClerk), intent(inout) :: this

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test kAlphaAnalogClerk for a 1 Cycle Batch case
  !!
@Test
  subroutine test1CycleBatch(this)
    class(test_kAlphaAnalogClerk), intent(inout) :: this
    type(scoreMemory)                            :: mem
    type(particleDungeon)                        :: pit
    type(particle)                               :: p
    real(defReal)                                :: k, STDk, alpha, STDa
    class(tallyResult), allocatable              :: res
    real(defReal),parameter :: TOL = 1.0E-9

    ! Initialise objects
    call pit % init(4)
    call mem % init(4_longInt,1)
    call this  % clerk % setMemAddress(1_longInt)

    ! Start cycle 1
    p % w = 1000.0_defReal
    call pit % detain(p)
    call this % clerk % reportCycleStart(pit, mem)

    ! End cycle 1
    call pit % release(p)
    p % w = 1200.0_defReal
    call pit % detain(p)
    pit % k_eff = ONE

    call this % clerk % reportCycleEnd(pit,mem)
    call mem % closeCycle(0.8_defReal)

    ! Start cycle 2
    call pit % release(p)
    p % w = 1000.0_defReal
    call pit % detain(p)
    call this % clerk % reportCycleStart(pit, mem)

    ! End cycle 2
    call pit % release(p)
    p % w = 900.0_defReal
    call pit % detain(p)
    pit % k_eff = 1.2_defReal

    call this % clerk % reportCycleEnd(pit,mem)
    call mem % closeCycle(0.8_defReal)

    ! Validate results

    ! Directly from memory
    call mem % getResult(k, STDk, 1_longInt)
    call mem % getResult(alpha, STDa, 2_longInt)
    @assertEqual(1.1400_defReal, k, TOL, '1 Cycle Batch, k from memory:')
    @assertEqual(0.0600_defReal, STDk, TOL, '1 Cycle Batch, k STD from memory:')
    @assertEqual(1.2480_defReal, alpha, TOL, '1 Cycle Batch, alpha from memory:')
    @assertEqual(0.0480_defReal, STDa, TOL, '1 Cycle Batch, alpha STD from memory:')

    ! From result
    call this % clerk % getResult(res, mem)

    ! Note these quantities are cycle-wise, not accumulated
    select type(res)
      type is(kAlphaResult)
        @assertEqual(1.0800_defReal, res % k, TOL, '1 Cycle Batch, k from result:')
        @assertEqual(1.2960_defReal, res % alpha, TOL, '1 Cycle Batch, alpha from result:')

      class default
        @assertTrue(.false.,'Result is not a kAlphaResult')

    end select

  end subroutine test1CycleBatch


  !!
  !! Test kAlphaAnalogClerk for a 2 Cycle Batch case
  !!
@Test
  subroutine test2CycleBatch(this)
    class(test_kAlphaAnalogClerk), intent(inout) :: this
    type(scoreMemory)                            :: mem
    type(particleDungeon)                        :: pit
    type(particle)                               :: p
    real(defReal)                                :: k, STDk, alpha, STDa
    class(tallyResult), allocatable              :: res
    real(defReal),parameter :: TOL = 1.0E-9

    ! Initialise objects
    call pit % init(4)
    call mem % init(2_longInt,1, batchSize = 2 )
    call this  % clerk % setMemAddress(1_longInt)

    ! Start cycle 1
    p % w = 500.0_defReal
    call pit % detain(p)
    call this % clerk % reportCycleStart(pit, mem)

    ! End cycle 1
    call pit % release(p)
    p % w = 500.0_defReal
    call pit % detain(p)
    pit % k_eff = ONE

    call this % clerk % reportCycleEnd(pit,mem)
    call mem % closeCycle(0.8_defReal)

    ! Start cycle 2
    call pit % release(p)
    p % w = 500.0_defReal
    call pit % detain(p)
    call this % clerk % reportCycleStart(pit, mem)

    ! End cycle 2
    call pit % release(p)
    p % w = 700.0_defReal
    call pit % detain(p)
    pit % k_eff = ONE

    call this % clerk % reportCycleEnd(pit,mem)
    call mem % closeCycle(0.8_defReal)

    ! Start cycle 3
    call pit % release(p)
    p % w = 500.0_defReal
    call pit % detain(p)
    call this % clerk % reportCycleStart(pit, mem)

    ! End cycle 3
    call pit % release(p)
    p % w = 400.0_defReal
    call pit % detain(p)
    pit % k_eff = 1.2_defReal

    call this % clerk % reportCycleEnd(pit,mem)
    call mem % closeCycle(0.8_defReal)

    ! Start cycle 4
    call pit % release(p)
    p % w = 500.0_defReal
    call pit % detain(p)
    call this % clerk % reportCycleStart(pit, mem)

    ! End cycle 4
    call pit % release(p)
    p % w = 500.0_defReal
    call pit % detain(p)
    pit % k_eff = 1.2_defReal

    call this % clerk % reportCycleEnd(pit,mem)
    call mem % closeCycle(0.8_defReal)


    ! Validate results

    ! Directly from memory
    call mem % getResult(k, STDk, 1_longInt)
    call mem % getResult(alpha, STDa, 2_longInt)
    @assertEqual(1.1400_defReal, k, TOL, '2 Cycle Batch, k from memory:')
    @assertEqual(0.0600_defReal, STDk, TOL, '2 Cycle Batch, k STD from memory:')
    @assertEqual(1.2480_defReal, alpha, TOL, '2 Cycle Batch, alpha from memory:')
    @assertEqual(0.0480_defReal, STDa, TOL, '2 Cycle Batch, alpha STD from memory:')

    ! From result
    call this % clerk % getResult(res, mem)

    select type(res)
      type is(kAlphaResult)
        @assertEqual(1.0800_defReal, res % k, TOL, '2 Cycle Batch, k from result:')
        @assertEqual(1.2960_defReal, res % alpha, TOL, '2 Cycle Batch, alpha STD from result:')

      class default
        @assertTrue(.false.,'Result is not a kAlphaResult')

    end select

  end subroutine test2CycleBatch

  !!
  !! Test functions: print & getSize
  !!
@Test
  subroutine testMisc(this)
    class(test_kAlphaAnalogClerk), intent(inout) :: this
    type(outputFile)                             :: out
    type(scoreMemory)                            :: mem

    ! Initialise objects
    call mem % init(2_longInt,1)
    call this  % clerk % setMemAddress(1_longInt)

    ! Test getting size
    @assertEqual(2, this % clerk % getSize(),'Test getSize() :')

    ! Test output printing correctness
    call out % init('dummyPrinter', fatalErrors = .false.)
    call this % clerk % print(out, mem)
    @assertTrue(out % isValid(), 'Test print():')

  end subroutine testMisc


end module kAlphaAnalogClerk_test
