module multiScatterMG_test

  use numPrecision
  use endfConstants
  use RNG_class,            only : RNG
  use dictionary_class,     only : dictionary
  use dictDeck_class,       only : dictDeck
  use reactionHandle_inter, only : reactionHandle
  use reactionMG_inter,     only : reactionMG, reactionMG_CptrCast
  use multiScatterMG_class, only : multiScatterMG, multiScatterMG_TptrCast, multiScatterMG_CptrCast
  use funit

  implicit none

  !! Test parameters
  real(defReal),dimension(4),parameter :: P0   = [1.3_defReal, 0.7_defReal, 0.3_defReal, 4.0_defReal ]
  real(defReal),dimension(4),parameter :: prod = [1.1_defReal, 1.05_defReal, ONE, ONE]


contains

  !!
  !! Tests pointer casting and basic (deterministic) functionality
  !! Does NOT test sampling procedure!
  !!
@Test
  subroutine multiplicative_scattering_test()
    type(multiScatterMG), target   :: reaction
    class(reactionHandle), pointer :: handlePtr
    class(reactionMG), pointer     :: mgPtr
    class(multiScatterMG), pointer :: multiCPtr
    type(multiScatterMG), pointer  :: multiTPtr
    type(dictionary),target       :: dictT
    type(dictDeck)                :: data
    type(RNG)                     :: rand
    real(defReal),parameter :: TOL = 1.0E-6_defReal

    ! Set pointers
    handlePtr => reaction
    mgPtr => null()
    multiCPtr => null()
    multiTPtr => null()

    ! Uncorrelated Reaction class cast
    mgPtr => reactionMG_CptrCast(reaction)
    @assertTrue(associated(mgPtr, reaction))

    ! Multiplicative Scattering class cast
    multiCPtr => multiScatterMG_CptrCast(reaction)
    @assertTrue(associated(multiCPtr, reaction))

    ! Multiplicative Scattering type cast
    multiTPtr => multiScatterMG_TptrCast(reaction)
    @assertTrue(associated(multiTPtr, reaction))

    ! Initialise multiplicative scattering
    call dictT % init(3)
    call dictT % store('numberOfGroups',2)
    call dictT % store('P0',P0)
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
    @assertEqual(ONE, reaction % production(2,1), TOL)
    @assertEqual(1.05_defReal, reaction % production(1,2), TOL)

    ! Clean memory
    call dictT % kill()
    call reaction % kill()

  end subroutine multiplicative_scattering_test


end module multiScatterMG_test
