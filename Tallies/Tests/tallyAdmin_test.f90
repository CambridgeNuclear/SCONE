module tallyAdmin_test

  use dictionary_class,      only : dictionary
  use dictParser_func,       only : charToDict
  use funit
  use keffAnalogClerk_class, only : keffResult
  use numPrecision
  use particle_class,        only : particle
  use particleDungeon_class, only : particleDungeon
  use tallyAdmin_class,      only : tallyAdmin
  use tallyResult_class,     only : tallyResult

  implicit none

  character(*), parameter :: TALLY_DEF = "clerk1 {type keffAnalogClerk;} &
                                         &norm clerk1; normVal 1.0; batchSize 2;"

contains

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test combining a normalisation clerk with batchSize > 1.
  !!
  !! Scores are folded into reduced score bins only on batch-closing cycles,
  !! hence the normalisation score must be read only on those cycles as well.
  !! Reading it every cycle finds empty bins in between and kills the run with 
  !! 'Normalisation score ... is 0' on the very first active cycle.
  !! The definition below follows the sample dictionary in the tallyAdmin
  !! documentation (normalisation clerk together with a batch size above 1).
  !!
  !! Two identical batches of two cycles with equal start and end populations
  !! must give a k-eff estimate of exactly 1.0 with zero spread.
  !!
@Test
  subroutine testBatchedNormalisation()
    class(tallyResult), allocatable :: res
    integer(shortInt)               :: i
    type(dictionary)                :: dict
    type(particle)                  :: p
    type(particleDungeon)           :: pop
    type(tallyAdmin)                :: tally
    real(defReal), parameter        :: TOL = 1.0E-9_defReal

    ! Build tally admin from definition.
    call charToDict(dict, TALLY_DEF)
    call tally % init(dict)
    call dict % kill()

    ! Populate the dungeon: 5 particles of weight 2 -> population weight 10.
    call pop % init(10)
    p % w = 2.0_defReal
    do i = 1, 5
      call pop % detain(p)

    end do

    ! Run 4 cycles = 2 batches with identical start and end populations.
    do i = 1, 4
      call tally % reportCycleStart(pop)
      call tally % reportCycleEnd(pop)

    end do

    ! k-eff estimate must be ONE with an STD of ZERO since the batches are identical.
    call tally % getResult(res, 'clerk1')

    select type(res)
      type is(keffResult)
        @assertEqual(ONE, res % keff(1), TOL, 'Batched normalisation, keff:')
        @assertEqual(ZERO, res % keff(2), TOL, 'Batched normalisation, keff STD:')

      class default
        @assertTrue(.false., 'Result is not a keffResult.')

    end select

    ! Clean up.
    call tally % kill()
    call pop % kill()

  end subroutine testBatchedNormalisation

end module tallyAdmin_test