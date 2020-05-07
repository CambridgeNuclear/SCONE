module pureCapture_test

  use numPrecision
  use RNG_class,                    only : RNG
  use reactionHandle_inter,         only : reactionHandle
  use uncorrelatedReactionCE_inter, only : uncorrelatedReactionCE, uncorrelatedReactionCE_CptrCast
  use pureCapture_class,            only : pureCapture, pureCapture_TptrCast
  use dictDeck_class,               only : dictDeck

  use pFUnit_mod

  implicit none


contains

  !!
  !! Integration test of elasticScattering reaction
  !! Tests:
  !!   -> building from ACE
  !!   -> pointer Casting
  !!   -> probability of Scattering
  !!
  !! Does NOT verify correctness of the sampling procedure
  !!
@Test
  subroutine testPureCaptureReaction()
    type(pureCapture), target             :: reaction
    class(reactionHandle),pointer         :: handlePtr
    class(uncorrelatedReactionCE),pointer :: unCorrPtr
    type(pureCapture),pointer             :: pureCapPtr
    type(dictDeck)                        :: fakeDeck
    type(RNG)                             :: rand
    real(defReal)                         :: E_out, mu, phi
    real(defReal),parameter :: TOL = 1.0E-6_defReal

    ! Set pointers
    handlePtr => reaction
    unCorrPtr => null()
    pureCapPtr => null()

    ! Uncorrelated Reaction cast
    unCorrPtr => uncorrelatedReactionCE_CptrCast(reaction)
    @assertTrue(associated(unCorrPtr, reaction))

    ! Elastic Scattering type cast
    pureCapPtr => pureCapture_TptrCast(reaction)
    @assertTrue(associated(pureCapPtr, reaction))

    ! Initialise -> not needed really
    call pureCapPtr % init(fakeDeck, 2)

    ! Test functionality
    ! Release & Delay rate sampling
    @assertEqual(ZERO, pureCapPtr % release(1.0_defReal))
    @assertEqual(ZERO, pureCapPtr % releasePrompt(2.0_defReal))
    @assertEqual(ZERO, pureCapPtr % releaseDelayed(2.0_defReal))
    @assertFalse(pureCapPtr % inCMframe())

    ! Verify probability distribution
    @assertEqual(ONE, pureCapPtr % probOf(ONE, ZERO, 1.0_defReal, 1.0_defReal))
    @assertEqual(ONE, pureCapPtr % probOf(ONE, TWO_PI, 1.0_defReal, 1.0_defReal))
    @assertEqual(ZERO, pureCapPtr % probOf(ONE, ZERO, 1.1_defReal, 1.0_defReal))
    @assertEqual(ZERO, pureCapPtr % probOf(ONE, 1.0_defReal, 1.0_defReal, 1.0_defReal))
    @assertEqual(ZERO, pureCapPtr % probOf(0.0_defReal, ZERO, 1.0_defReal, 1.0_defReal))

    ! Sample outgoing
    call pureCapPtr % sampleOut(mu, phi, E_out, 8.0_defReal, rand)
    @assertEqual(ONE, mu, TOL)
    @assertEqual(ZERO, phi, TOL)
    @assertEqual(8.0_defReal, E_out, TOL)

    ! Clean
    call reaction % kill()

  end subroutine testPureCaptureReaction


end module pureCapture_test
