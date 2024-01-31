module keffAnalogClerk_test

  use numPrecision
  use tallyResult_class,     only : tallyResult
  use keffAnalogClerk_class, only : keffAnalogClerk, keffResult
  use particle_class,        only : particle
  use particleDungeon_class, only : particleDungeon
  use dictionary_class,      only : dictionary
  use scoreMemory_class,     only : scoreMemory
  use outputFile_class,      only : outputFile
  use funit

  implicit none

@testCase
  type, extends(TestCase) :: test_keffAnalogClerk
    private
    type(keffAnalogClerk) :: clerk
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_keffAnalogClerk

contains

  !!
  !! Sets up test_keffAnalogClerk object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_keffAnalogClerk), intent(inout) :: this
    type(dictionary)                           :: dict
    character(nameLen)                         :: name

    call dict % init(2)
    call this % clerk % init(dict, name)

    call dict % kill()

  end subroutine setUp

  !!
  !! Kills test_keffAnalogClerk object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_keffAnalogClerk), intent(inout) :: this

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test keffAnalogClerk for a 1 Cycle Batch case
  !!
@Test
  subroutine test1CycleBatch(this)
    class(test_keffAnalogClerk), intent(inout) :: this
    type(scoreMemory)                          :: mem
    type(particleDungeon)                      :: pit
    type(particle)                             :: p
    real(defReal)                              :: k, STD
    class(tallyResult), allocatable            :: res
    real(defReal),parameter :: TOL = 1.0E-9

    ! Initialise objects
    call pit % init(4)
    call mem % init(2_longInt,1)
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

    call mem % reduceBins()
    call this % clerk % reportCycleEnd(pit, mem)
    call mem % closeCycle(0.8_defReal)

    ! Validate results

    ! Directly from memory
    call mem % getResult(k, STD, 1_longInt)
    @assertEqual(1.1400_defReal, k, TOL, '1 Cycle Batch, keff from memory:')
    @assertEqual(0.0600_defReal, STD, TOL, '1 Cycle Batch, keff STD from memory:')

    ! From result
    call this % clerk % getResult(res, mem)

    select type(res)
      type is(keffResult)
        @assertEqual(1.1400_defReal, res % keff(1), TOL, '1 Cycle Batch, keff from result:')
        @assertEqual(0.0600_defReal, res % keff(2), TOL, '1 Cycle Batch, keff STD from result:')

      class default
        @assertTrue(.false.,'Result is not a keffResult')

    end select

  end subroutine test1CycleBatch


  !!
  !! Test keffAnalogClerk for a 2 Cycle Batch case
  !!
@Test
  subroutine test2CycleBatch(this)
    class(test_keffAnalogClerk), intent(inout) :: this
    type(scoreMemory)                          :: mem
    type(particleDungeon)                      :: pit
    type(particle)                             :: p
    real(defReal)                              :: k, STD
    class(tallyResult), allocatable            :: res
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
    call mem % getResult(k, STD, 1_longInt)
    @assertEqual(1.1400_defReal, k, TOL, '1 Cycle Batch, keff from memory:')
    @assertEqual(0.0600_defReal, STD, TOL, '1 Cycle Batch, keff STD from memory:')

    ! From result
    call this % clerk % getResult(res, mem)

    select type(res)
      type is(keffResult)
        @assertEqual(1.1400_defReal, res % keff(1), TOL, '1 Cycle Batch, keff from result:')
        @assertEqual(0.0600_defReal, res % keff(2), TOL, '1 Cycle Batch, keff STD from result:')

      class default
        @assertTrue(.false.,'Result is not a keffResult')

    end select

  end subroutine test2CycleBatch

  !!
  !! Test functions: print & getSize
  !!
@Test
  subroutine testMisc(this)
    class(test_keffAnalogClerk), intent(inout) :: this
    type(outputFile)                           :: out
    type(scoreMemory)                          :: mem

    ! Initialise objects
    call mem % init(2_longInt,1)
    call this  % clerk % setMemAddress(1_longInt)

    ! Test getting size
    @assertEqual(1, this % clerk % getSize(),'Test getSize() :')

    ! Test output printing correctness
    call out % init('dummyPrinter', fatalErrors = .false.)
    call this % clerk % print(out, mem)
    @assertTrue(out % isValid(), 'Test print():')

  end subroutine testMisc


end module keffAnalogClerk_test
