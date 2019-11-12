module neutronScattering_iTest

  use numPrecision
  use endfConstants
  use RNG_class,                    only : RNG
  use reactionHandle_inter,         only : reactionHandle
  use uncorrelatedReactionCE_inter, only : uncorrelatedReactionCE, uncorrelatedReactionCE_CptrCast
  use neutronScatter_class,         only : neutronScatter, neutronScatter_TptrCast
  use aceCard_class,                only : aceCard

  use pFUnit_mod
  implicit none


contains

  !!
  !! Integration test of elasticScattering reaction
  !! Tests:
  !!   -> building from ACE
  !!   -> pointer Casting
  !!
  !! There were some problems to establish same normalisation wrt JANIS
  !! Leave probability test for now!
  !!
  !! Does NOT verify correctness of the sampling procedure
  !!
@Test
  subroutine testNeutronScatteringReaction()
    type(neutronScatter), target          :: reaction
    class(reactionHandle),pointer         :: handlePtr
    class(uncorrelatedReactionCE),pointer :: unCorrPtr
    type(neutronScatter),pointer          :: scatterPtr
    type(aceCard)                         :: ACE
    type(RNG)                             :: rand
    real(defReal),parameter :: TOL = 1.0E-6_defReal

    ! Set pointers
    handlePtr => reaction
    unCorrPtr => null()
    scatterPtr => null()

    ! Uncorrelated Reaction cast
    unCorrPtr => uncorrelatedReactionCE_CptrCast(reaction)
    @assertTrue(associated(unCorrPtr, reaction))

    ! Elastic Scattering type cast
    scatterPtr => neutronScatter_TptrCast(reaction)
    @assertTrue(associated(scatterPtr, reaction))

    ! Build ACE library
    call ACE % readFromFile('./IntegrationTestFiles/8016JEF311.ace', 1)

    ! Build reaction object
    call reaction % init(ACE, N_2N)

    ! Verify simple functionality
    @assertTrue(reaction % inCMFrame())
    @assertEqual(TWO,  reaction % release(7.0_defReal))
    @assertEqual(TWO,  reaction % releasePrompt(1.0E-6_defReal))
    @assertEqual(ZERO, reaction % releaseDelayed(1.5E-3_defReal))

    ! Clean
    call reaction % kill()

  end subroutine testNeutronScatteringReaction


end module neutronScattering_iTest
