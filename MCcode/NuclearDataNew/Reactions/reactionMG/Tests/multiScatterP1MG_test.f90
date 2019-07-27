module multiScatterP1MG_test

  use numPrecision
  use endfConstants
  use RNG_class,              only : RNG
  use dictionary_class,       only : dictionary
  use dictDeck_class,         only : dictDeck
  use reactionHandle_inter,   only : reactionHandle
  use reactionMG_inter,       only : reactionMG, reactionMG_CptrCast
  use multiScatterP1MG_class, only : multiScatterP1MG, multiScatterP1MG_TptrCast
  use pfUnit_mod

  implicit none

  !! Test parameters
  real(defReal),dimension(4),parameter :: P0   = [1.3_defReal, 0.7_defReal, 0.3_defReal, 4.0_defReal ]
  real(defReal),dimension(4),parameter :: P1   = [0.5_defReal, 0.1_defReal, ZERO, ZERO ]
  real(defReal),dimension(4),parameter :: prod = [1.1_defReal, 1.05_defReal, ONE, ONE]

contains

  !!
  !! Tests pointer casting and basic (deterministic) functionality
  !! Does NOT test sampling procedure!
  !!
  !! NOTE: Mostly copy of test for multiScatterMG. Significant functionality is
  !!   shared, but it is not necessarly so. Thus test is duplicated.
  !!
@Test
  subroutine multiplicative_p1_scattering_test()
    type(multiScatterP1MG), target   :: reaction
    class(reactionHandle), pointer   :: handlePtr
    class(reactionMG), pointer       :: mgPtr
    type(multiScatterP1MG), pointer  :: multiTPtr
    type(dictionary),target          :: dictT
    type(dictDeck)                   :: data
    type(RNG)                        :: rand
    real(defReal),parameter :: TOL = 1.0E-6_defReal

    ! Set pointers
    handlePtr => reaction
    mgPtr => null()
    multiTPtr => null()

    ! Uncorrelated Reaction class cast
    mgPtr => reactionMG_CptrCast(reaction)
    @assertTrue(associated(mgPtr, reaction))

    ! Multiplicative Scattering type cast
    multiTPtr => multiScatterP1MG_TptrCast(reaction)
    @assertTrue(associated(multiTPtr, reaction))

    ! Initialise multiplicative scattering
    call dictT % init(4)
    call dictT % store('numberOfGroups',2)
    call dictT % store('P0',P0)
    call dictT % store('P1',P1)
    call dictT % store('scatteringMultiplicity', prod)
    data % dict => dictT
    call reaction % init(data, anyScatter)

    ! Test misc functionality
    @assertEqual(ZERO, reaction % releaseDelayed(1), TOL)
    @assertEqual(ZERO, reaction % sampleDelayRate(1, rand), TOL)

    ! Test average release
    @assertEqual(ZERO, reaction % releasePrompt(-1), TOL)
    @assertEqual(ZERO, reaction % release(-1), TOL)
    @assertEqual(ZERO, reaction % releasePrompt(172), TOL)
    @assertEqual(ZERO, reaction % release(42), TOL)
    @assertEqual(ONE, reaction % release(2), TOL)
    @assertEqual(1.0825_defReal, reaction % releasePrompt(1), TOL)

    ! Test getting total scattering Xss
    @assertEqual(2.0_defReal, reaction % scatterXS(1), TOL)
    @assertEqual(4.3_defReal, reaction % scatterXS(2), TOL)

    ! Test getting Group-To-Group production
    @assertEqual(1.1_defReal, reaction % production(1,1), TOL)
    @assertEqual(1.05_defReal, reaction % production(2,1), TOL)
    @assertEqual(ONE, reaction % production (1,2), TOL)

    ! Test that normalisation of P1 coefficients is OK
    @assertEqual(ZERO, reaction % P1(2,2), TOL)
    @assertEqual(ZERO, reaction % P1(1,2), TOL)
    @assertEqual(0.3846153846_defReal, reaction % P1(1,1), TOL)
    @assertEqual(0.1428571429_defReal, reaction % P1(2,1), TOL)

    ! Clean memory
    call dictT % kill()
    call reaction % kill()

  end subroutine multiplicative_p1_scattering_test
    
end module multiScatterP1MG_test
