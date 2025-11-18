module particle_test
  use numPrecision
  use universalVariables
  use particle_class, only : particle, particleState, P_NEUTRON, P_PHOTON, P_PRECURSOR, verifyType
  use funit

  implicit none

@TestCase
  type, extends(TestCase) :: test_particle
    type(particle) :: p_MG
    type(particle) :: p_CE
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_particle

  ! Test paramethers
  real(defReal),dimension(3),parameter :: r0 = [-1.3_defReal, 97.7_defReal, -9.0_defReal]
  real(defReal),dimension(3),parameter :: u0 = [1.0_defReal, 0.0_defReal, 0.0_defReal]
  real(defReal),parameter              :: t0 = 13.0_defReal
  real(defReal),parameter              :: w0 = 0.7_defReal
  real(defReal),parameter              :: E0 = 0.00034_defReal
  integer(shortInt),parameter          :: G0 = 3

  real(defReal),dimension(3),parameter :: lev1_offset = [0.3_defReal, -13.5_defReal, 1.0_defReal]
  real(defReal),dimension(3),parameter :: lev2_offset = [-3.4_defReal, 3.0_defReal, -8.0_defReal]

  integer(shortInt),parameter          :: lev1_uni     = 5
  integer(shortInt),parameter          :: lev1_uniRoot = 67
  integer(shortInt),parameter          :: lev2_uni     = 3
  integer(shortInt),parameter          :: lev2_uniRoot = 32

contains

  !!
  !! Set-up particle for testing
  !!
  subroutine setUp(this)
    class(test_particle), intent(inout) :: this

    ! Build MG particle
    call this % p_MG % build(r0, u0, G0, w0, t0)
    this % p_MG % type = P_NEUTRON

    ! Build CE particle
    call this % p_CE % build(r0, u0, E0, w0)
    this % p_CE % type = P_PHOTON

    ! Add 2 levels to the particles
    call this % p_MG % coords % addLevel()
    this % p_MG % coords % lvl(2) % r = r0 - lev1_offset
    this % p_MG % coords % lvl(2) % dir = u0

    call this % p_CE % coords % addLevel()
    this % p_CE % coords % lvl(2) % r = r0 - lev1_offset
    this % p_CE % coords % lvl(2) % dir = u0

    call this % p_MG % coords % addLevel()
    this % p_MG % coords % lvl(3) % r = r0 - lev1_offset - lev2_offset
    this % p_MG % coords % lvl(3) % dir = u0

    call this % p_CE % coords % addLevel()
    this % p_CE % coords % lvl(3) % r = r0 - lev1_offset - lev2_offset
    this % p_CE % coords % lvl(3) % dir = u0

    ! Set MatIdx
    this % p_MG % coords % matIdx = 7
    this % p_CE % coords % matIdx = 7

    ! Set uniqueID
    this % p_MG % coords % uniqueID = 34
    this % p_CE % coords % uniqueID = 34

    ! NOTE: THIS MAY BREAK AT SOME POINT
    ! HAND SET SOME COORD PARAMETERS
    ! BREAKS ENCAPSULATION (It is broken anyway)

    ! Level 1
    this % p_MG % coords % lvl(1) % uniIdx    = 1
    this % p_MG % coords % lvl(1) % uniRootID = 0
    this % p_MG % coords % lvl(1) % localID   = 1
    this % p_MG % coords % lvl(1) % cellIdx   = 3

    this % p_CE % coords % lvl(1) % uniIdx    = 1
    this % p_CE % coords % lvl(1) % uniRootID = 0
    this % p_CE % coords % lvl(1) % localID   = 1
    this % p_CE % coords % lvl(1) % cellIdx   = 3

    ! Level 2
    this % p_MG % coords % lvl(2) % uniIdx    = lev1_uni
    this % p_MG % coords % lvl(2) % uniRootID = lev1_uniRoot
    this % p_MG % coords % lvl(2) % localID   = 4
    this % p_MG % coords % lvl(2) % cellIdx   = 2

    this % p_CE % coords % lvl(2) % uniIdx    = lev1_uni
    this % p_CE % coords % lvl(2) % uniRootID = lev1_uniRoot
    this % p_CE % coords % lvl(2) % localID   = 4
    this % p_CE % coords % lvl(2) % cellIdx   = 2

    ! Level 3
    this % p_MG % coords % lvl(3) % uniIdx    = lev2_uni
    this % p_MG % coords % lvl(3) % uniRootID = lev2_uniRoot
    this % p_MG % coords % lvl(3) % localID   = 2
    this % p_MG % coords % lvl(3) % cellIdx   = 8

    this % p_CE % coords % lvl(3) % uniIdx    = lev2_uni
    this % p_CE % coords % lvl(3) % uniRootID = lev2_uniRoot
    this % p_CE % coords % lvl(3) % localID   = 2
    this % p_CE % coords % lvl(3) % cellIdx   = 8

  end subroutine setUp

  !!
  !! Kill particle after test
  !!
  subroutine tearDown(this)
    class(test_particle), intent(inout) :: this
    ! DO NOTHING
  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Tests begin here
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test initialisation correctness
  !!
@test
  subroutine correctInitialisation(this)
    class(test_particle), intent(inout) :: this

    ! Test energies
    @assertEqual(E0, this % p_CE % E, 'CE energy value')
    @assertEqual(G0, this % p_MG % G, 'MG group number')

    ! Test time
    @assertEqual(ZERO, this % p_CE % time ,'Default time initialisation')
    @assertEqual(t0, this % p_MG % time ,'Time initialisation')

    ! Test weight
    @assertEqual(w0, this % p_CE % w, 'CE weight initialisation')
    @assertEqual(w0, this % p_CE % w0, 'CE initial weigth initialisation')

    @assertEqual(w0, this % p_MG % w, 'MG weight initialisation')
    @assertEqual(w0, this % p_MG % w0, 'MG initial weigth initialisation')

    ! Test isMG flag
    @assertFalse( this % p_CE % isMG, 'isMG flag for CE particle')
    @assertTrue( this % p_MG % isMG, 'isMG flag for MG particle')

    ! Test isDead flag
    @assertFalse( this % p_CE % isDead, 'isDead flag CE particle')
    @assertFalse( this % p_MG % isDead, 'isDead flag MG particle')

    ! Test timeMax default value
    @assertEqual(-INF, this % p_CE % timeMax, 'timeMax initialises to -INF by default')
    @assertEqual(-INF, this % p_MG % timeMax, 'timeMax initialises to -INF by default')
    
    ! Test lambda default value
    @assertEqual(INF, this % p_CE % lambda, 'lambda initialises to INF by default')
    @assertEqual(INF, this % p_MG % lambda, 'lambda initialises to INF by default')

  end subroutine correctInitialisation

  !!
  !! Test access to positions
  !!
@test
  subroutine testPositionAccess(this)
    class(test_particle), intent(inout) :: this
    real(defReal),dimension(3)          :: r_global
    real(defReal),dimension(3)          :: r_lev1
    real(defReal),dimension(3)          :: r_lev2
    real(defReal),dimension(3)          :: r_lev3

    ! Get position at all levels
    r_global = this % p_CE % rGlobal()

    r_lev1 = this % p_CE % rLocal(1)
    r_lev2 = this % p_CE % rLocal(2)
    r_lev3 = this % p_CE % rLocal()

    ! Verify correctness
    @assertEqual(r0,r_global,'Global Position Wrong')
    @assertTrue(r_global == r_lev1, 'Level 1 r is not equal global r')
    @assertEqual(r0-lev1_offset, r_lev2, 'Level 2 position is wrong')
    @assertEqual(r0-lev1_offset-lev2_offset, r_lev3, 'Level 3 position is wrong')

  end subroutine testPositionAccess

  !!
  !! Test access to directions
  !!
@test
  subroutine testDirectionAccess(this)
    class(test_particle), intent(inout) :: this
    real(defReal),dimension(3)          :: u_global
    real(defReal),dimension(3)          :: u_lev1
    real(defReal),dimension(3)          :: u_lev2
    real(defReal),dimension(3)          :: u_lev3

    ! Get direction at all levels
    u_global = this % p_CE % dirGlobal()

    u_lev1 = this % p_CE % dirLocal(1)
    u_lev2 = this % p_CE % dirLocal(2)
    u_lev3 = this % p_CE % dirLocal()

    ! Verify correctness
    @assertEqual(u0,u_global,'Global Direction Wrong')
    @assertTrue(u_global == u_lev1, 'Level 1 dir is not equal global dir')
    @assertEqual(u0, u_lev2, 'Level 2 Direction is wrong')
    @assertEqual(u0, u_lev3, 'Level 3 Direction is wrong')

  end subroutine testDirectionAccess

  !!
  !! Tests access procedures to extra parameters. TODO: add uniqueID
  !!
@test
  subroutine testMiscAccess(this)
    class(test_particle), intent(inout) :: this
    integer(shortInt)                   :: matIdx, cellIdx, uniIdx

    ! Verify nesting level
    @assertEqual(3, this % p_CE % nesting(), 'Nesting Level')

    ! Level 3
    matIdx  = this % p_CE % matIdx()
    cellIdx = this % p_CE % getCellIdx()
    uniIdx  = this % p_CE % getUniIdx()

    ! Verify correctness
    @assertEqual(7, matIdx, 'Material Index')
    @assertEqual(8, cellIdx, 'Cell Index. Deepest level.')
    @assertEqual(3, uniIdx, 'Universe Index. Deepest level.')

    ! Level 2
    cellIdx = this % p_CE % getCellIdx(2)
    uniIdx  = this % p_CE % getUniIdx(2)

    ! Verify correctness
    @assertEqual(2, cellIdx, 'Cell Index. Level 2.')
    @assertEqual(5, uniIdx, 'Universe Index. Level 2.')

    ! Level 1
    cellIdx = this % p_CE % getCellIdx(1)
    uniIdx  = this % p_CE % getUniIdx(1)

    ! Verify correctness
    @assertEqual(3, cellIdx, 'Cell Index. Level 1.')
    @assertEqual(1, uniIdx, 'Universe Index. Level 1.')

    ! Verify getting speed
    @assertEqual(lightSpeed, this % p_CE % getSpeed())

  end subroutine testMiscAccess

  !!
  !! Test material setting outside geometry
  !!
@test
  subroutine testSetMatIdx(this)
    class(test_particle), intent(inout) :: this

    call this % p_CE % setMatIdx(3)

    @assertEqual(3, this % p_CE % matIdx())

  end subroutine testSetMatIdx

  !!
  !! Test movement procedures
  !!
@test
  subroutine testMovementProcedures(this)
    class(test_particle), intent(inout) :: this
    real(defReal)                       :: dist = 1.0_defReal
    real(defReal),dimension(3)          :: r0_lvl1 ,r0_lvl2, r0_lvl3
    real(defReal),dimension(3)          :: r_lvl1, r_lvl2, r_lvl3

    ! Move local on lowest level
    call this % p_CE % moveLocal(dist,3)

    r_lvl1 = this % p_CE % rLocal(1)
    r_lvl2 = this % p_CE % rLocal(2)
    r_lvl3 = this % p_CE % rLocal()

    ! Calculate expected position
    r0_lvl1 = r0 + u0 * dist
    r0_lvl2 = r0 - lev1_offset + u0 * dist
    r0_lvl3 = r0 - lev1_offset - lev2_offset + u0 * dist

    ! Verify position
    @assertEqual(r0_lvl1, r_lvl1, 'Level 1 position. Local Movement.')
    @assertEqual(r0_lvl2, r_lvl2, 'Level 2 position. Local Movement.')
    @assertEqual(r0_lvl3, r_lvl3, 'Level 3 position. Local Movement.')
    @assertTrue(this % p_CE % coords % isPlaced(), 'Particle is placed in the geometry')

    ! Move on global level
    call this % p_CE % moveGlobal(dist)
    r_lvl1  = this % p_CE % rGlobal()
    r0_lvl1 = r0_lvl1 + u0 * dist

    @assertEqual(r0_lvl1, r_lvl1, 'Global position after global movement')
    @assertTrue(this % p_CE % coords % isAbove(), 'Particle is above geometry')

    ! Teleport on global level
    r0_lvl1 = [-1.0_defReal, 6.0_defReal, -14.6868_defReal]
    call this % p_CE % teleport(r0_lvl1)
    r_lvl1 = this % p_CE % rGlobal()

    @assertEqual(r0_lvl1, r_lvl1, 'Global position after global teleport')
    @assertTrue(this % p_CE % coords % isAbove(), 'Particle is above geometry')

  end subroutine testMovementProcedures

  !!
  !!
  !!
@test
  subroutine testRotationProcedures(this)
    class(test_particle), intent(inout) :: this
    real(defReal), dimension(3)         :: dir
    real(defReal)                       :: tol

    dir = [4.0_defReal, -5.0_defReal, 1.0_defReal]
    dir = dir / norm2(dir)

    ! Test point procedure
    call this % p_CE % point(dir)

    @assertEqual(dir, this % p_CE % dirGlobal(), 'Global direction after pointing')
    @assertEqual(dir, this % p_CE % dirLocal(1),' Local direction after pointing. Level 1')
    @assertEqual(dir, this % p_CE % dirLocal(2),' Local direction after pointing. Level 2')
    @assertEqual(dir, this % p_CE % dirLocal(),' Local direction after pointing. Level 3')
    @assertTrue(this % p_CE % coords % isPlaced(), 'Particle is still placed in geometry')

    ! Test rotation
    call this % p_CE % rotate(0.3_defReal, 1.3_defReal)

    ! Rotation was performed with an independent, verified MATLAB implementation
    dir = [0.927517049363521_defReal, 0.312003100719400_defReal,   -0.205830484334726_defReal]
    tol = 50.0 * epsilon(dir)
    @assertEqual(dir, this % p_CE % dirGlobal(), tol, 'Global direction after rotation')
    @assertEqual(dir, this % p_CE % dirLocal(1), tol,'Local direction after rotation. Level 1')
    @assertEqual(dir, this % p_CE % dirLocal(2), tol,'Local direction after rotation. Level 2')
    @assertEqual(dir, this % p_CE % dirLocal(), tol,'Local direction after rotation. Level 3')
    @assertTrue(this % p_CE % coords % isPlaced(),'Particle is still placed in geometry')

    ! Test taking above geometry
    call this % p_CE % takeAboveGeom()
    @assertTrue(this % p_CE % coords % isAbove(), 'Particle is above geometry')

  end subroutine testRotationProcedures

  !!
  !! Test that assignments between particle and particleState are correct
  !!
@test
  subroutine testParticleStateAssignments(this)
    class(test_particle), intent(inout) :: this
    type(particleState)                 :: mg_coord
    type(particleState)                 :: ce_coord
    type(particle)                      :: mg_p
    type(particle)                      :: ce_p

    ! Copy particle to phaseCoords
    mg_coord = this % p_MG
    ce_coord = this % p_CE

    ! Copy phaseCoords to particle
    mg_p = mg_coord
    ce_p = ce_coord

    ! Verify correctness (Use only MG execept where it differs from CE: energy and isMG flag)

    ! Compare global position and direction
    @assertEqual(r0, mg_coord % r, 'Global position. PhaseCoord')
    @assertEqual(r0, mg_p % rGlobal(), 'Global position. Particle from PhaseCoord')

    @assertEqual(u0, mg_coord % dir, 'Global direction, PhaseCoord')
    @assertEqual(u0, mg_p % dirGLobal(), 'Global direction. Particle from PhaseCoord')

    ! Compare weight
    @assertEqual(w0, mg_coord % wgt, 'Weight. PhaseCoord')
    @assertEqual(w0, mg_p % w, ' Weight. Particle from PhaseCoord')
    @assertEqual(w0, mg_p % w0, 'Starting weight. Particle from PhaseCoord')

    ! Compare time
    @assertEqual(t0, mg_coord % time, 'Time. PhaseCoord')
    @assertEqual(t0, mg_p % time, 'Time. Particle from PhaseCoord')

    ! Compare energy
    @assertEqual(G0, mg_coord % G, 'Energy group. PhaseCoord')
    @assertEqual(E0, ce_coord % E, 'Energy. PhaseCoord')
    @assertTrue(mg_coord % isMG, 'isMG flag for MG. PhaseCoord')
    @assertFalse(ce_coord % isMG, 'isMG flag for CE. PhaseCoord')

    @assertEqual(G0, mg_p % G, 'Energy group. Particle from PhaseCoord')
    @assertEqual(E0, ce_p % E, 'Energy. Particle from PhaseCoord')
    @assertTrue(mg_p % isMG, 'isMG flag for MG. Particle from PhaseCoord')
    @assertFalse(ce_p % isMG, 'isMG flag for CE. Particle from PhaseCoord')

    ! Compare Type
    @assertEqual(P_PHOTON, ce_p % type, 'Type of particle has changed')
    @assertEqual(P_NEUTRON, mg_p % type, 'Type of particle has changed')

    ! Verify material, cell IDXs and unique ID
    @assertEqual(7, mg_coord % matIdx, 'Material index')
    @assertEqual(8, mg_coord% cellIdx, 'Cell index')
    @assertEqual(34, mg_coord % uniqueID, 'Unique ID')
  end subroutine testParticleStateAssignments

!  !!
!  !! Test that assignment between particle and particleState is correct
!  !!   DOES NOT RETEST ASSIGNMENT OF phaseCoords parts
!  !!
!@test
!  subroutine testParticleStateAssignment(this)
!    class(test_particle), intent(inout) :: this
!    type(particleState)                 :: state
!
!    state = this % p_MG
!
!    ! Verify
!
!
!  end subroutine testParticleStateAssignment

  !!
  !! Test state saving procedures
  !!
@test
  subroutine testStateSaving(this)
    class(test_particle), intent(inout) :: this
    type(particleState)                 :: stateRef
    logical(defBool)                    :: isCorrect

    ! **** Do CE particle first*****************************
    ! Save reference
    stateRef = this % p_CE

    ! Check prehistory saving
    call this % p_CE % savePreHistory()
    isCorrect = this % p_CE % preHistory == stateRef
    @assertTrue(isCorrect,'preHistory check CE')

    ! Check pretransition saving
    call this % p_CE % savePreTransition()
    isCorrect = this % p_CE % preTransition == stateRef
    @assertTrue(isCorrect,'preTransition check CE')

    ! Check prePath saving
    call this % p_CE % savePrePath()
    isCorrect = this % p_CE % prePath == stateRef
    @assertTrue(isCorrect,'prePath check CE')

    ! Check precollision saving
    call this % p_CE % savePreCollision()
    isCorrect = this % p_CE % preCollision == stateRef
    @assertTrue(isCorrect,'preCollision check CE')

    ! **** Do MG particle **********************************
    ! Save reference
    stateRef = this % p_MG

    ! Check prehistory saving
    call this % p_MG % savePreHistory()
    isCorrect = this % p_MG % preHistory == stateRef
    @assertTrue(isCorrect,'preHistory check MG')

    ! Check pretransition saving
    call this % p_MG % savePreTransition()
    isCorrect = this % p_MG % preTransition == stateRef
    @assertTrue(isCorrect,'preTransition check MG')

    ! Check prePath saving
    call this % p_MG % savePrePath()
    isCorrect = this % p_MG % prePath == stateRef
    @assertTrue(isCorrect,'prePath check MG')

    ! Check precollision saving
    call this % p_MG % savePreCollision()
    isCorrect = this % p_MG % preCollision == stateRef
    @assertTrue(isCorrect,'preCollision check MG')

  end subroutine testStateSaving

  !!
  !! Test type verification
  !!
@Test
  subroutine testParticleTypeVerification(this)
    class(test_particle), intent(inout) :: this

    @assertTrue( verifyType(P_NEUTRON), 'Particle Neutron')
    @assertTrue( verifyType(P_PHOTON), 'Particle Photon')
    @assertTrue( verifyType(P_PRECURSOR), 'Particle Precursor')
    @assertFalse( verifyType(-876864), 'Invalid particle type parameter')

  end subroutine testParticleTypeVerification

  !!
  !! Test precursor procedures
  !!
  subroutine testPrecursorProcedures(this)
    class(test_particle), intent(inout) :: this
    type(particle)                      :: that
    real(defReal), parameter            :: TOL = 1.0E-6

    ! Precursor type check
    this % p_CE % type = P_NEUTRON
    @assertFalse( this % p_CE % isPrecursor(), 'Particle Neutron')
    this % p_CE % type = P_PRECURSOR
    @assertTrue( this % p_CE % isPrecursor(), 'Particle Precursor')

    ! Convert precursor to neutron
    this % p_CE % lambda = ONE
    call this % p_CE % emitDelayedNeutron()
    @assertFalse( this % p_CE % isPrecursor(), 'Particle Neutron')
    @assertEqual( this % p_CE % lambda, INF)

    ! Test forced precursor decay
    this % p_CE % lambda = ONE
    this % p_CE % type = P_PRECURSOR
    this % p_CE % w = ONE
    this % p_CE % time = ONE
    call this % p_CE % forcedPrecursorDecay(TWO, TWO, that)
    @assertTrue( this % p_CE % isPrecursor(), 'Particle Precursor')
    @assertFalse( that % isPrecursor(), 'Particle Neutron')
    @assertEqual( that % w, 0.735758882_defReal, TOL)

  end subroutine testPrecursorProcedures

end module particle_test
