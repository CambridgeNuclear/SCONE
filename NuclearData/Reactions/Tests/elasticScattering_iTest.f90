module elasticScattering_iTest

  use numPrecision
  use endfConstants
  use RNG_class,                    only : RNG
  use reactionHandle_inter,         only : reactionHandle
  use uncorrelatedReactionCE_inter, only : uncorrelatedReactionCE, uncorrelatedReactionCE_CptrCast
  use elasticNeutronScatter_class,  only : elasticNeutronScatter, elasticNeutronScatter_TptrCast
  use aceCard_class,                only : aceCard
  use funit
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
    unCorrPtr => uncorrelatedReactionCE_CptrCast(reaction)
    @assertTrue(associated(unCorrPtr, reaction))

    ! Elastic Scattering type cast
    elasticScatterPtr => elasticNeutronScatter_TptrCast(reaction)
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

    ! Varify probability distributions
    ! Isotropic range
    @assertEqual(ONE/TWO/TWO_PI, reaction % probOf(0.5_defReal, 2.0_defReal, 0.6E-6_defReal, 0.6E-6_defReal), TOL)

    ! High energy
    ! It is based on regression really but verified against JANIS polynomial representation
    ! with error of about 1.0E-3 absolute, which is accaptable
    ! It was immposible to access polynomial representation directly, and there might be additional
    ! Issues due to tabular representation of the data.
    @assertEqual(0.4273399E-01_defReal, reaction % probOf(0.7_defReal, 2.0_defReal, 0.36_defReal, 0.36_defReal), TOL)
    @assertEqual(0.1593347_defReal, reaction % probOf(0.75911_defReal, 2.0_defReal, 3.94_defReal, 3.94_defReal), TOL)
    @assertEqual(0.7943654E-01_defReal, reaction % probOf(-0.3_defReal, 2.0_defReal, 8.9201_defReal, 8.9201_defReal), TOL)

    ! Test invalid angle ranges
    @assertEqual(ZERO, reaction % probOf(1.1_defReal, 2.0_defReal, ONE, ONE), TOL)
    @assertEqual(ZERO, reaction % probOf(0.7_defReal, -2.0_defReal, ONE, ONE), TOL)
    @assertEqual(ZERO, reaction % probOf(0.7_defReal, 2.0_defReal, -ONE, -ONE), TOL)

    ! Clean
    call reaction % kill()

  end subroutine testElasticNeutronScatteringReaction

  !!
  !! Test elasticScattering for a data with LOCB == 0
  !! Use Ta-126 from JEFF 3.1.1
  !!
  !! Does not test sampling
  !!
@Test
  subroutine testElasticNeutronScattering_isotropic()
    type(elasticNeutronScatter)  :: reaction
    type(aceCard)                :: ACE
    real(defReal),parameter :: TOL = 1.0E-6_defReal

    ! Build ACE library
    call ACE % readFromFile('./IntegrationTestFiles/52126JEF311.ace', 1)

    ! Build reaction object
    call reaction % init(ACE, N_N_ELASTIC)

    ! Test probability distribution
    @assertEqual(ONE/TWO/TWO_PI, reaction % probOf(0.5_defReal, 2.0_defReal, 6.7_defReal, 6.7_defReal), TOL)

    ! Test invalid angle aranges
    @assertEqual(ZERO, reaction % probOf(1.1_defReal, 2.0_defReal, ONE, ONE), TOL)
    @assertEqual(ZERO, reaction % probOf(0.7_defReal, -2.0_defReal, ONE, ONE), TOL)
    @assertEqual(ZERO, reaction % probOf(0.7_defReal, 2.0_defReal, -ONE, -ONE), TOL)

  end subroutine testElasticNeutronScattering_isotropic

end module elasticScattering_iTest
