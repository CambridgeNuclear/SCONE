module imcMaterialSource_class

  use numPrecision
  use endfConstants
  use universalVariables
  use genericProcedures,       only : fatalError, rotateVector
  use dictionary_class,        only : dictionary
  use RNG_class,               only : RNG

  use particle_class,          only : particle, particleState, P_PHOTON, P_MATERIAL
  use particleDungeon_class,   only : particleDungeon
  use source_inter,            only : source, kill_super => kill

  use geometry_inter,          only : geometry
  use geometryGrid_class,      only : geometryGrid
  use IMCMaterial_inter,       only : IMCMaterial, IMCMaterial_CptrCast
  use nuclearDataReg_mod,      only : ndReg_getIMCMG => getIMCMG
  use nuclearDatabase_inter,   only : nuclearDatabase
  use mgIMCDatabase_inter,     only : mgIMCDatabase
  use materialMenu_mod,        only : mm_nMat    => nMat, &
                                      mm_matName => matName

  use simulationTime_class

  implicit none
  private

  ! Calculation type
  integer(shortInt), parameter :: IMC = 1, ISMC = 2

  !!
  !! Source for uniform generation of photons within a material for IMC and ISMC calculations
  !!
  !! Angular distribution is isotropic
  !!
  !! Private members:
  !!   isMG     -> is the source multi-group? (default = .true.)
  !!   bottom   -> Bottom corner (x_min, y_min, z_min)
  !!   top      -> Top corner (x_max, y_max, z_max)
  !!   G        -> Energy group
  !!   pType    -> P_PHOTON for IMC, P_MATERIAL for ISMC
  !!   bounds   -> Bounds of geometry
  !!   calcType -> IMC or ISMC, changes type of material to be sampled
  !!
  !! Interface:
  !!   source_inter Interface
  !!
  !! SAMPLE INPUT:
  !!   matSource { type imcMaterialSource; calcType IMC; }
  !!
  type, public,extends(source) :: imcMaterialSource
    private
    logical(defBool)                :: isMG     = .true.
    real(defReal), dimension(3)     :: bottom   = ZERO
    real(defReal), dimension(3)     :: top      = ZERO
    integer(shortInt)               :: G        = 0
    integer(shortInt)               :: pType    = P_PHOTON
    real(defReal), dimension(6)     :: bounds   = ZERO
    integer(shortInt)               :: calcType = IMC
  contains
    procedure :: init
    procedure :: append
    procedure :: sampleParticle
    procedure, private :: sampleIMC
    procedure :: kill
  end type imcMaterialSource

contains

  !!
  !! Initialise material Source
  !!
  !! See source_inter for details
  !!
  subroutine init(self, dict, geom)
    class(imcMaterialSource), intent(inout)     :: self
    class(dictionary), intent(in)            :: dict
    class(geometry), pointer, intent(in)     :: geom
    character(100), parameter :: Here = 'init (imcMaterialSource_class.f90)'

    call dict % getOrDefault(self % G, 'G', 1)

    ! Provide geometry info to source
    self % geom => geom

    ! Set bounding region
    self % bounds = self % geom % bounds()

    ! Select calculation type - Automatically added to dict in implicitPhysicsPackage
    call dict % getOrDefault(self % calcType, 'calcType', IMC)
    select case(self % calcType)
      case(IMC)
        self % pType = P_PHOTON
      case(ISMC)
        self % pType = P_MATERIAL
      case default
        call fatalError(Here, 'Unrecognised calculation type. Should be "IMC" or "ISMC"')
    end select

  end subroutine init


  !!
  !! Generate N particles to add to a particleDungeon without overriding particles already present.
  !! Note that energy here refers to energy weight rather than group.
  !!
  !! Args:
  !!   dungeon [inout] -> particle dungeon to be added to
  !!   n [in]          -> number of particles to place in dungeon
  !!   rand [inout]    -> particle RNG object
  !!
  !! Result:
  !!   A dungeon populated with N particles sampled from the source, plus particles
  !!   already present in dungeon
  !!
  subroutine append(self, dungeon, N, rand)
    class(imcMaterialSource), intent(inout)    :: self
    type(particleDungeon), intent(inout)    :: dungeon
    integer(shortInt), intent(in)           :: N
    class(RNG), intent(inout)               :: rand
    real(defReal), dimension(6)             :: bounds
    integer(shortInt)                       :: matIdx, i, Ntemp, G
    real(defReal)                           :: energy, totalEnergy
    type(RNG)                               :: pRand
    class(mgIMCDatabase), pointer           :: nucData
    class(geometry), pointer                :: geom
    character(100), parameter               :: Here = "append (imcMaterialSource_class.f90)"

    ! Get pointer to appropriate nuclear database
    nucData => ndReg_getIMCMG()
    if(.not.associated(nucData)) call fatalError(Here, 'Failed to retrieve Nuclear Database')

    ! Obtain total energy
    if (self % calcType == IMC) then
      totalEnergy = nucData % getEmittedRad()
    else
      totalEnergy = nucData % getMaterialEnergy()
    end if

    ! Loop through materials
    do matIdx = 1, mm_nMat()

      ! Get energy to be emitted from material matIdx
      if (self % calcType == IMC) then
        energy = nucData % getEmittedRad(matIdx)
      else
        energy = nucData % getMaterialEnergy(matIdx)
      end if

      ! Choose particle numbers in proportion to material energy
      if (energy > ZERO) then
        Ntemp = int(N * energy / totalEnergy)
        ! Enforce at least 1 particle
        if (Ntemp == 0) Ntemp = 1

        ! Set bounds for sampling
        geom => self % geom
        select type(geom)
          class is(geometryGrid)
            bounds = geom % matBounds(matIdx)
          class default
            bounds = self % bounds
          end select

        ! Find energy per particle
        energy = energy / Ntemp

        ! Sample particles
        !$omp parallel
        pRand = rand
        !$omp do private(pRand, G)
        do i=1, Ntemp
          call pRand % stride(i)
          G = nucData % sampleEnergyGroup(matIdx, pRand)
          call dungeon % detain(self % sampleIMC(pRand, matIdx, energy, G, bounds))
        end do
        !$omp end do
        !$omp end parallel

      end if
    end do
   
  end subroutine append

  !!
  !! Should not be called
  !!
  function sampleParticle(self, rand) result(p)
    class(imcMaterialSource), intent(inout) :: self
    class(RNG), intent(inout)            :: rand
    type(particleState)                  :: p
    character(100), parameter :: Here = 'sampleParticle (imcMaterialSource_class.f90)'

    ! Should not be called, useful to have extra inputs so use sampleIMC instead
    call fatalError(Here, 'Should not be called, sampleIMC should be used instead.')

    ! Avoid compiler warning
    p % G = self % G

  end function sampleParticle

  !!
  !! Sample particle's phase space co-ordinates
  !!
  !! Args:
  !!   rand [in]   -> RNG
  !!   matIdx [in] -> index of material being sampled from
  !!   energy [in] -> energy-weight of sampled particle
  !!   G [in]      -> energy group of sampled particle
  !!   bounds [in] -> bounds for position search, will be bounds of entire geometry if using
  !!                  geometryStd, and bounds of single material if using geometryGrid
  !!
  function sampleIMC(self, rand, targetMatIdx, energy, G, bounds) result(p)
    class(imcMaterialSource), intent(inout)    :: self
    class(RNG), intent(inout)               :: rand
    integer(shortInt), intent(in)           :: targetMatIdx
    real(defReal), intent(in)               :: energy
    integer(shortInt), intent(in)           :: G
    real(defReal), dimension(6), intent(in) :: bounds
    type(particleState)                     :: p
    real(defReal), dimension(3)             :: bottom, top, r, dir, rand3
    real(defReal)                           :: mu, phi
    integer(shortInt)                       :: i, matIdx, uniqueID
    character(100), parameter :: Here = 'sampleIMC (imcMaterialSource_class.f90)'

    ! Sample particle position
    bottom = bounds(1:3)
    top    = bounds(4:6)
    i = 0
    rejection:do
      rand3(1) = rand % get()
      rand3(2) = rand % get()
      rand3(3) = rand % get()
      r = (top - bottom) * rand3 + bottom

      ! Find material under position
      call self % geom % whatIsAt(matIdx, uniqueID, r)

      ! Exit if in desired material
      if (matIdx == targetMatIdx) exit rejection

      ! Protect against infinite loop
      i = i+1
      if (i > 10000) call fatalError(Here, '10,000 failed attempts in rejection sampling')

    end do rejection

    ! Sample direction - chosen uniformly inside unit sphere
    mu = 2 * rand % get() - 1
    phi = rand % get() * 2*pi
    dir(1) = mu
    dir(2) = sqrt(1-mu**2) * cos(phi)
    dir(3) = sqrt(1-mu**2) * sin(phi)

    ! Sample time uniformly within time step
    p % time = time % stepStart + timeStep() * rand % get()

    ! Assign basic phase-space coordinates
    p % matIdx   = matIdx
    p % uniqueID = uniqueID
    p % r        = r
    p % dir      = dir
    p % G        = G
    p % isMG     = .true.
    p % wgt      = energy
    p % type     = self % pType

  end function sampleIMC

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(imcMaterialSource), intent(inout) :: self

    call kill_super(self)

    self % isMG   = .true.
    self % bounds = ZERO
    self % G      = 0

  end subroutine kill

end module imcMaterialSource_class
