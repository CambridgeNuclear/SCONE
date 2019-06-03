module elasticScattering_iTest

  use numPrecision
  use endfConstants
  use RNG_class,                    only : RNG
  use reactionHandle_inter,         only : reactionHandle
  use uncorrelatedReactionCE_inter, only : uncorrelatedReactionCE, uncorrelatedReactionCE_ptrCast
  use elasticNeutronScatter_class,  only : elasticNeutronScatter, elasticNeutronScatter_ptrCast
  use aceCard_class,                only : aceCard
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
  subroutine testElasticNeutronScatteringReaction()
    type(elasticNeutronScatter), target   :: reaction
    class(reactionHandle),pointer         :: handlePtr
    class(uncorrelatedReactionCE),pointer :: unCorrPtr
    type(elasticNeutronScatter),pointer   :: elasticScatterPtr
    type(aceCard)                         :: ACE
    type(RNG)                             :: rand
    real(defReal),parameter :: TOL = 1.0E-6_defReal

    ! Set pointers
    handlePtr => reaction
    unCorrPtr => null()
    elasticScatterPtr => null()

    ! Uncorrelated Reaction cast
    unCorrPtr => uncorrelatedReactionCE_ptrCast(reaction)
    @assertTrue(associated(unCorrPtr, reaction))

    ! Elastic Scattering type cast
    elasticScatterPtr => elasticNeutronScatter_ptrCast(reaction)
    @assertTrue(associated(elasticScatterPtr, reaction))

    ! Build ACE library
    call ACE % readFromFile('./IntegrationTestFiles/8016JEF311.ace', 1)

    ! Build reaction object
    call reaction % init(ACE, N_N_ELASTIC)

    ! Verify simple functionality
    @assertTrue(reaction % inCMFrame())
    @assertEqual(ONE,  reaction % release(7.0_defReal))
    @assertEqual(ONE,  reaction % releasePrompt(1.0E-6_defReal))
    @assertEqual(ZERO, reaction % releaseDelayed(1.5E-3_defReal))
    @assertEqual(ZERO, reaction % sampleDelayRate(1.5E-3_defReal, rand))

    ! Varify probability distributions
    ! Isotropic range
    @assertEqual(ONE/TWO/TWO_PI, reaction % probOf(0.5_defReal, 2.0_defReal, 0.6E-6_defReal, 0.6E-6_defReal), TOL)

    ! High energy
    ! It is based on regression really but verified against JANIS polynomial representation
    ! with error of about 1.0E-3 absolute, which is accaptable
    ! It was immposible to access polynomial representation directly, and there might be additional
    ! Issues due to tabular representation of the data.
    @assertEqual(0.4257917E-01_defReal, reaction % probOf(0.7_defReal, 2.0_defReal, 0.36_defReal, 0.36_defReal), TOL)
    @assertEqual(0.1593337_defReal, reaction % probOf(0.75911_defReal, 2.0_defReal, 3.94_defReal, 3.94_defReal), TOL)
    @assertEqual(0.7934733E-01_defReal, reaction % probOf(-0.3_defReal, 2.0_defReal, 8.9201_defReal, 8.9201_defReal), TOL)

    ! Test invalid angle aranges
    @assertEqual(ZERO, reaction % probOf(1.1_defReal, 2.0_defReal, ONE, ONE), TOL)
    @assertEqual(ZERO, reaction % probOf(0.7_defReal, -2.0_defReal, ONE, ONE), TOL)
    @assertEqual(ZERO, reaction % probOf(0.7_defReal, 2.0_defReal, -ONE, -ONE), TOL)

    ! Clean
    call reaction % kill()

  end subroutine testElasticNeutronScatteringReaction


    
end module elasticScattering_iTest
