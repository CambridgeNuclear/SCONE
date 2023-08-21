module fissionMG_test

  use numPrecision
  use endfConstants
  use RNG_class,            only : RNG
  use dictionary_class,     only : dictionary
  use dictDeck_class,       only : dictDeck
  use reactionHandle_inter, only : reactionHandle
  use reactionMG_inter,     only : reactionMG, reactionMG_CptrCast
  use fissionMG_class,      only : fissionMG, fissionMG_TptrCast
  use funit

  implicit none

  !!
  !! Test data
  !!
  real(defReal),dimension(*),parameter :: nu = [2.3_defReal, 2.0_defReal, 1.3_defReal]
  real(defReal),dimension(*),parameter :: chi = [0.333333_defReal, 0.333333_defReal, 0.333334_defReal]



contains

  !!
  !! Tests pointer casting and basic (deterministic) functionality
  !! Does NOT test sampling procedure!
  !!
@Test
  subroutine fissionMG_Build_And_Functionality()
    type(fissionMG), target       :: reaction
    class(reactionHandle),pointer :: handlePtr
    class(reactionMG),pointer     :: mgPtr
    type(fissionMG),pointer       :: fissPtr
    type(dictionary),target       :: dictT
    type(dictionary),pointer      :: dictPtr
    type(dictDeck)                :: data
    type(RNG)                     :: rand
    real(defReal),parameter :: TOL = 1.0E-6_defReal

    ! Set pointers
    handlePtr => reaction
    mgPtr => null()
    fissPtr => null()

    ! Uncorrelated Reaction cast
    mgPtr => reactionMG_CptrCast(reaction)
    @assertTrue(associated(mgPtr, reaction))

    ! Elastic Scattering type cast
    fissPtr => fissionMG_TptrCast(reaction)
    @assertTrue(associated(fissPtr, reaction))

    ! Build dictionary for input
    call dictT % init(3)
    call dictT % store('numberOfGroups',3)
    call dictT % store('chi', chi)
    call dictT % store('nu',nu)

    ! Build data Deck and initialise
    data % dict => dictT
    call reaction % init(data, macroFission)

    ! Test Misc functionality
    @assertEqual(ZERO, reaction % releaseDelayed(1), TOL)
    @assertEqual(ZERO, reaction % sampleDelayRate(2, rand),TOL )

    ! Test Release
    @assertEqual(ZERO, reaction % releasePrompt(-2), TOL)
    @assertEqual(ZERO, reaction % release(6), TOL)
    @assertEqual(2.0_defReal, reaction % release(2), TOL)
    @assertEqual(1.3_defReal,reaction % releasePrompt(3), TOL)

    ! Clean
    call dictT % kill()
    call reaction % kill()

  end subroutine fissionMG_Build_And_Functionality
    
end module fissionMG_test
