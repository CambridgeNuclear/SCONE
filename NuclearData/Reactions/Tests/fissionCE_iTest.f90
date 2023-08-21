module fissionCE_iTest

  use numPrecision
  use endfConstants
  use RNG_class,                    only : RNG
  use reactionHandle_inter,         only : reactionHandle
  use uncorrelatedReactionCE_inter, only : uncorrelatedReactionCE, uncorrelatedReactionCE_CptrCast
  use fissionCE_class,              only : fissionCE, fissionCE_TptrCast
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
  subroutine testFissionCE()
    type(fissionCE), target               :: reaction
    class(reactionHandle),pointer         :: handlePtr
    class(uncorrelatedReactionCE),pointer :: unCorrPtr
    type(fissionCE),pointer               :: fissionPtr
    type(aceCard)                         :: ACE
    type(RNG)                             :: rand
    real(defReal),parameter :: TOL = 1.0E-6_defReal

    ! Set pointers
    handlePtr => reaction
    unCorrPtr => null()
    fissionPtr => null()

    ! Uncorrelated Reaction cast
    unCorrPtr => uncorrelatedReactionCE_CptrCast(reaction)
    @assertTrue(associated(unCorrPtr, reaction))

    ! Elastic Scattering type cast
    fissionPtr => fissionCE_TptrCast(reaction)
    @assertTrue(associated(fissionPtr, reaction))

    ! Build ACE library
    call ACE % readFromFile('./IntegrationTestFiles/92233JEF311.ace', 1)

    ! Build reaction object
    call reaction % init(ACE, N_FISSION)

    ! Test trivial functionality
    @assertFalse(reaction % inCMframe())

    ! Test neutron release
    @assertEqual(2.65431_defReal, reaction % release(1.6_defReal), TOL)
    @assertEqual(5.147534_defReal, reaction % release(17.0_defReal), TOL)
    @assertEqual(2.48098_defReal, reaction % releasePrompt(0.6E-6_defReal), TOL)
    @assertEqual(6.73E-3_defReal, reaction % releaseDelayed(0.6E-6_defReal), TOL)
    @assertEqual(0.0041725_defReal, reaction % releaseDelayed(17.0_defReal), TOL)

    ! Test probability density
    !@assertEqual(0.1618843E-01_defReal, reaction % probOf(0.7_defReal, 2.0_defReal, 0.1404_defReal, 2.0_defReal), TOL)
    !@assertEqual(0.1586432E-01, reaction % probOf(0.7_defReal, 2.0_defReal, 2.48077_defReal, 14.0_defReal), TOL)

    ! Clean
    call reaction % kill()

  end subroutine testFissionCE

end module fissionCE_iTest
