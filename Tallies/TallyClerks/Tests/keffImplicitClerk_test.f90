module keffImplicitClerk_test

  use numPrecision
  use tallyCodes
  use endfConstants
  use tallyResult_class,       only : tallyResult
  use keffAnalogClerk_class,   only : keffResult
  use keffImplicitClerk_class, only : keffImplicitClerk
  use particle_class,          only : particle
  use particleDungeon_class,   only : particleDungeon
  use dictionary_class,        only : dictionary
  use scoreMemory_class,       only : scoreMemory
  use outputFile_class,        only : outputFile
  use funit

  use testNeutronDatabase_class, only : testNeutronDatabase

  implicit none

@testCase
  type, extends(TestCase) :: test_keffImplicitClerk
    private
    type(keffImplicitClerk)   :: clerk
    type(testNeutronDatabase) :: nucData
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_keffImplicitClerk

contains

  !!
  !! Sets up test_keffImplicitClerk object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_keffImplicitClerk), intent(inout) :: this
    type(dictionary)                             :: dict
    character(nameLen)                           :: name

    call dict % init(2)
    call this % clerk % init(dict, name)
    call dict % kill()

    call this % nucData % build(ONE, captureXS = 2.0_defReal, fissionXS = ONE, nuFissionXS = 3.0_defReal)

  end subroutine setUp

  !!
  !! Kills test_keffImplicitClerk object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_keffImplicitClerk), intent(inout) :: this

    call this % nucData % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test for 1 cycle batches
  !!
@Test
  subroutine test1CycleBatch(this)
    class(test_keffImplicitClerk), intent(inout) :: this
    type(particle)                               :: p
    type(particleDungeon)                        :: pit
    type(scoreMemory)                            :: mem
    class(tallyResult),allocatable               :: res
    real(defReal), parameter                     :: TOL = 1.0E-9_defReal

    ! Configure memory
    call mem % init(10_longInt, 1)
    call this % clerk % setMemAddress(1_longInt)
    call pit % init(4)


    ! Configure particle
    p % fate = leak_FATE

    !*** Start cycle 1
    ! Score implicit reaction rates
    p % w = 0.7_defReal
    call this % clerk % reportInColl(p, this % nucData, mem, .false.)

    ! Score analog production
    p % preCollision % wgt = 0.1_defReal
    call this % clerk % reportOutColl(p, N_2N, 0.5_defReal, this % nucData, mem)

    ! Score leakage
    p % w = 0.3_defReal
    call this % clerk % reportHist(p, this % nucData, mem)

    ! End cycle
    call mem % reduceBins()
    call pit % detain(p)
    call this % clerk % closeCycle(pit, mem)
    call pit % release(p)
    call mem % closeCycle(ONE)

    !*** Start cycle 2
    ! Score implicit reaction rates
    p % w = 0.6_defReal
    call this % clerk % reportInColl(p, this % nucData, mem, .false.)

    ! Score analog production
    p % preCollision % wgt = 0.1_defReal
    call this % clerk % reportOutColl(p, N_2N, 0.5_defReal, this % nucData, mem)

    ! Score leakage
    p % w = 0.3_defReal
    call this % clerk % reportHist(p, this % nucData, mem)

    ! End cycle
    call mem % reduceBins()
    call pit % detain(p)
    call this % clerk % closeCycle(pit, mem)
    call pit % release(p)
    call mem % closeCycle(ONE)

    ! Verify result
    call this % clerk % getResult(res, mem)
    select type(res)
      type is(keffResult)
        @assertEqual(0.906521739130435_defReal, res % keff(1), TOL, '1 Cycle Batch, keff from result:')
        @assertEqual(0.006521739130435_defReal, res % keff(2), TOL, '1 Cycle Batch, keff STD from result:')

      class default
        @assertTrue(.false.,'Result is not a keffResult')

    end select

  end subroutine test1CycleBatch

  !!
  !! Test getSize() and print
  !!
@Test
  subroutine testMisc(this)
    class(test_keffImplicitClerk), intent(inout) :: this
    type(scoreMemory)                            :: mem
    type(outputFile)                             :: out

    ! Configure memory
    call mem % init(10_longInt, 1)
    call this % clerk % setMemAddress(1_longInt)
    call out % init('dummyPrinter', fatalErrors = .false.)

    ! Test getting size
    @assertEqual(5, this % clerk % getSize(),'Test getSize():')

    ! Test correctness of output calls
    call this % clerk % print(out, mem)
    @assertTrue(out % isValid(), 'Test print():')

  end subroutine testMisc

end module keffImplicitClerk_test
