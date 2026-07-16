module sixFactorClerk_test

  use numPrecision
  use tallyCodes
  use endfConstants
  use tallyResult_class,     only : tallyResult
  use particle_class,        only : particle
  use particleDungeon_class, only : particleDungeon
  use dictionary_class,      only : dictionary
  use scoreMemory_class,     only : scoreMemory
  use outputFile_class,      only : outputFile
  use funit
  use sixFactorClerk_class,  only : sixFactorClerk, FAST_FISSION_FACTOR, &
                                    FAST_NON_LEAKAGE, RESONANCE_ESCAPE_PROB, &
                                    THERMAL_NON_LEAKAGE, THERMAL_UTILISATION, &
                                    ETA_THERMAL, K_EFFECTIVE

  use testNeutronDatabase_class, only : testNeutronDatabase

  implicit none

@testCase
  type, extends(TestCase) :: test_sixFactorClerk
    private
    type(sixFactorClerk)      :: clerk
    type(testNeutronDatabase) :: nucData
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_sixFactorClerk

contains

  !!
  !! Sets up test_sixFactorClerk object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_sixFactorClerk), intent(inout) :: this
    type(dictionary)                          :: dict
    character(nameLen)                        :: name

    call dict % init(2)
    call this % clerk % init(dict, name)
    call dict % kill()

    call this % nucData % build(ONE, captureXS = 2.0_defReal, fissionXS = ONE, nuFissionXS = 3.0_defReal)

  end subroutine setUp

  !!
  !! Kills test_sixFactorClerk object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_sixFactorClerk), intent(inout) :: this

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
    class(test_sixFactorClerk), intent(inout) :: this
    type(particle)                            :: p
    type(particleDungeon)                     :: pit
    type(scoreMemory)                         :: mem
    class(tallyResult),allocatable            :: res
    real(defReal)                             :: epsilon, STDe, fastPNL, STDfPNL, resP, STDp, &
                                                 thermalPNL, STDtPNL, f, STDf, eta, STDeta, &
                                                 keff, STDk
    real(defReal), parameter                  :: TOL = 1.0E-6_defReal

    ! Configure memory
    call mem % init(20_longInt, 1)
    call this % clerk % setMemAddress(1_longInt)
    call pit % init(4)

    p % isMG = .false.

    !*** Start cycle 1
    ! Score implicit fast reaction rates
    p % w = 0.7_defReal
    p % E  = 0.5_defReal
    call this % clerk % reportInColl(p, this % nucData, mem, .false.)

    ! Score analog production
    p % preCollision % wgt = 0.1_defReal
    call this % clerk % reportOutColl(p, N_2N, 0.5_defReal, this % nucData, mem)
    
    ! Score fast leakage
    p % w = 0.8_defReal
    p % fate = leak_FATE
    call this % clerk % reportHist(p, this % nucData, mem)
    
    ! Score implicit thermal reaction rates
    p % w = 0.6_defReal
    p % E = 1.0E-7_defReal
    call this % clerk % reportInColl(p, this % nucData, mem, .false.)

    ! Score thermal leakage
    p % w = 0.3_defReal
    p % fate = leak_FATE
    call this % clerk % reportHist(p, this % nucData, mem)

    ! End cycle
    call mem % reduceBins()
    call pit % detain(p)
    call this % clerk % closeCycle(pit, mem)
    call pit % release(p)
    call mem % closeCycle(ONE)
    
    !*** Start cycle 2
    ! Score implicit fast reaction rates
    p % isMG = .false.
    p % w = 1.2_defReal
    p % E  = 5.0_defReal
    call this % clerk % reportInColl(p, this % nucData, mem, .false.)

    ! Score analog production
    p % preCollision % wgt = 0.2_defReal
    call this % clerk % reportOutColl(p, N_2N, 0.5_defReal, this % nucData, mem)
    
    ! Score fast leakage
    p % w = 0.5_defReal
    p % fate = leak_FATE
    call this % clerk % reportHist(p, this % nucData, mem)
    
    ! Score implicit thermal reaction rates
    p % w = 0.2_defReal
    p % E = 1.0E-8_defReal
    call this % clerk % reportInColl(p, this % nucData, mem, .false.)

    ! Score thermal leakage
    p % w = 0.99_defReal
    p % fate = leak_FATE
    call this % clerk % reportHist(p, this % nucData, mem)
    
    ! End cycle
    call mem % reduceBins()
    call pit % detain(p)
    call this % clerk % closeCycle(pit, mem)
    call pit % release(p)
    call mem % closeCycle(ONE)
    
    ! Get results
    call mem % getResult(epsilon, STDe, this % clerk % getMemAddress() + FAST_FISSION_FACTOR)
    call mem % getResult(fastPNL, STDfPNL, this % clerk % getMemAddress() + FAST_NON_LEAKAGE)
    call mem % getResult(resP, STDp, this % clerk % getMemAddress() + RESONANCE_ESCAPE_PROB)
    call mem % getResult(thermalPNL, STDtPNL, this % clerk % getMemAddress() + THERMAL_NON_LEAKAGE)
    call mem % getResult(f, STDf, this % clerk % getMemAddress() + THERMAL_UTILISATION)
    call mem % getResult(eta, STDeta, this % clerk % getMemAddress() + ETA_THERMAL)
    call mem % getResult(keff, STDk, this % clerk % getMemAddress() + K_EFFECTIVE)

    ! Verify result
    @assertEqual(4.77777778_defReal, epsilon, TOL, '1 Cycle Batch, epsilon from memory:')
    @assertEqual(0.876063269_defReal, fastPNL, TOL, '1 Cycle Batch, fast PNL from memory:')
    @assertEqual(0.403179191_defReal, resP, TOL, '1 Cycle Batch, resonance escape probability from memory:')
    @assertEqual(0.617250674_defReal, thermalPNL, TOL, '1 Cycle Batch, thermal PNL from memory:')
    @assertEqual(1.0_defReal, f, TOL, '1 Cycle Batch, thermal utilisation from memory:')
    @assertEqual(1.0_defReal, eta, TOL, '1 Cycle Batch, eta from memory:')
    @assertEqual(0.786643234_defReal, keff, TOL, '1 Cycle Batch, keff from memory:')
    @assertEqual(2.555555556_defReal, STDe, TOL, '1 Cycle Batch, epsilon std from memory:')
    @assertEqual(0.036063269_defReal, STDfPNL, TOL, '1 Cycle Batch, fast PNL std from memory:')
    @assertEqual(0.096820809_defReal, STDp, TOL, '1 Cycle Batch, resonance escape probability std from memory:')
    @assertEqual(0.239892183_defReal, STDtPNL, TOL, '1 Cycle Batch, thermal PNL std from memory:')
    @assertEqual(0.0_defReal, STDf, TOL, '1 Cycle Batch, thermal utilisation std from memory:')
    @assertEqual(0.0_defReal, STDeta, TOL, '1 Cycle Batch, eta std from memory:')
    @assertEqual(0.013356766_defReal, STDk, TOL, '1 Cycle Batch, keff std from memory:')

  end subroutine test1CycleBatch

  !!
  !! Test getSize() and print
  !!
@Test
  subroutine testMisc(this)
    class(test_sixFactorClerk), intent(inout) :: this
    type(scoreMemory)                         :: mem
    type(outputFile)                          :: out

    ! Configure memory
    call mem % init(10_longInt, 1)
    call this % clerk % setMemAddress(1_longInt)
    call out % init('dummyPrinter', fatalErrors = .false.)

    ! Test getting size
    @assertEqual(14, this % clerk % getSize(),'Test getSize():')

    ! Test correctness of output calls
    call this % clerk % print(out, mem)
    @assertTrue(out % isValid(), 'Test print():')

  end subroutine testMisc

end module sixFactorClerk_test
