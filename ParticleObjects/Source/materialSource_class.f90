module materialSource_class

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
  use IMCMaterial_inter,       only : IMCMaterial, IMCMaterial_CptrCast
  use nuclearDataReg_mod,      only : ndReg_getIMCMG => getIMCMG
  use nuclearDatabase_inter,   only : nuclearDatabase
  use mgIMCDatabase_inter,     only : mgIMCDatabase
  use materialMenu_mod,        only : mm_nMat    => nMat, &
                                      mm_matName => matName

  use simulationTime_class

  implicit none
  private

  ! Position sampling method
  integer(shortInt), parameter :: REJ = 1, FAST = 2
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
  !!   latPitch -> Pitch of lattice (if using a lattice geom)
  !!   latSizeN -> Lattice dimensions (if using a lattice geom)
  !!   G        -> Energy group
  !!   pType    -> P_PHOTON for IMC, P_MATERIAL for ISMC
  !!   bounds   -> Bounds of geometry
  !!   method   -> REJ uses rejection sampling for position (VERY slow for many materials)
  !!            -> FAST samples only within bounds of each material so no rejection needed,
  !!               currently only works for lattics geometry (hence lattice settings above)
  !!   calcType -> IMC or ISMC, changes type of material to be sampled
  !!
  !! Interface:
  !!   source_inter Interface
  !!
  !! SAMPLE INPUT:
  !!   matSource { type materialSource; method fast; }
  !!
  type, public,extends(source) :: materialSource
    private
    logical(defBool)                :: isMG     = .true.
    real(defReal), dimension(3)     :: bottom   = ZERO
    real(defReal), dimension(3)     :: top      = ZERO
    real(defReal), dimension(3)     :: latPitch = ZERO
    integer(shortInt), dimension(3) :: latSizeN = 0
    integer(shortInt)               :: G        = 0
    integer(shortInt)               :: pType    = P_PHOTON
    real(defReal), dimension(6)     :: bounds   = ZERO
    integer(shortInt)               :: method   = REJ
    integer(shortInt)               :: calcType = IMC
  contains
    procedure :: init
    procedure :: append
    procedure :: sampleParticle
    procedure, private :: sampleIMC
    procedure, private :: getMatBounds
    procedure :: kill
  end type materialSource

contains

  !!
  !! Initialise material Source
  !!
  !! See source_inter for details
  !!
  subroutine init(self, dict, geom)
    class(materialSource), intent(inout)     :: self
    class(dictionary), intent(in)            :: dict
    class(geometry), pointer, intent(in)     :: geom
    character(nameLen)                       :: method
    character(100), parameter :: Here = 'init (materialSource_class.f90)'

    call dict % getOrDefault(self % G, 'G', 1)

    ! Provide geometry info to source
    self % geom => geom

    ! Set bounding region
    self % bounds = self % geom % bounds()

    ! Select method for position sampling
    call dict % getOrDefault(method, 'method', 'rejection')
    select case(method)
      case('rejection')
        self % method = REJ

      case('fast')
        self % method = FAST
        ! Get lattice dimensions
        self % latSizeN = self % geom % latSizeN()
        self % latPitch = (self % bounds(4:6) - self % bounds(1:3)) / self % latSizeN

      case default
        call fatalError(Here, 'Unrecognised method. Should be "rejection" or "fast"')
    end select

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
  !! Generate N particles to add to a particleDungeon without overriding
  !! particles already present.
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
    class(materialSource), intent(inout)    :: self
    type(particleDungeon), intent(inout)    :: dungeon
    integer(shortInt), intent(in)           :: N
    class(RNG), intent(inout)               :: rand
    real(defReal), dimension(6)             :: bounds
    integer(shortInt)                       :: matIdx, i, Ntemp, G
    real(defReal)                           :: energy, totalEnergy
    type(RNG)                               :: pRand
    class(mgIMCDatabase), pointer           :: nucData
    character(100), parameter               :: Here = "append (materialSource_class.f90)"

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
        if (self % method == FAST) then
          bounds = self % getMatBounds(matIdx)
        else
          bounds = self % bounds
        end if

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
    class(materialSource), intent(inout) :: self
    class(RNG), intent(inout)            :: rand
    type(particleState)                  :: p
    character(100), parameter :: Here = 'sampleParticle (materialSource_class.f90)'

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
  !!                  rejection sampling method, and bounds of single material if using fast
  !!
  function sampleIMC(self, rand, targetMatIdx, energy, G, bounds) result(p)
    class(materialSource), intent(inout)    :: self
    class(RNG), intent(inout)               :: rand
    integer(shortInt), intent(in)           :: targetMatIdx
    real(defReal), intent(in)               :: energy
    integer(shortInt), intent(in)           :: G
    real(defReal), dimension(6), intent(in) :: bounds
    type(particleState)                     :: p
    real(defReal), dimension(3)             :: bottom, top, r, dir, rand3
    real(defReal)                           :: mu, phi
    integer(shortInt)                       :: i, matIdx, uniqueID
    character(100), parameter :: Here = 'sampleIMC (materialSource_class.f90)'

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

      ! Should exit immediately if using fast method as bounds should contain only matIdx
      if (self % method == FAST) call fatalError(Here, 'Fast sourcing returned incorrect material')

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
  !! Get location of material in lattice for position sampling
  !!
  !! Note that this may be incorrect depending on how lattice input is given, this function
  !! assumes that geometry has been generated by discretiseGeom_class.f90
  !!
  !! Args:
  !!   matIdx [in]     -> matIdx for which to calculate bounds
  !!   matBounds [out] -> boundary of lattice cell, [xmin,ymin,zmin,xmax,ymax,zmax]
  !!
  !! TODO:
  !!   Would be nice to have most of this in a geometry module
  !!
  function getMatBounds(self, matIdx) result(matBounds)
    class(materialSource), intent(inout)     :: self
    integer(shortInt), intent(in)            :: matIdx
    real(defReal), dimension(6)              :: matBounds
    integer(shortInt), dimension(3)          :: ijk
    integer(shortInt)                        :: i, latIdFlipped
    character(nameLen)                       :: matName
    character(100), parameter                :: Here = 'getMatBounds (materialSourceClass.f90)'

    ! Extract lattice position from mat name (e.g. "m106 -> 106")
    ! This is different from localID in latUniverse_class as is counting from a different
    ! corner (see get_ijk function description below)
    matName = mm_matName(matIdx)
    read (matName(2:), '(I10)') latIdFlipped

    ! Set bounds of lattice cell containing matIdx
    ijk = get_ijk(latIdFlipped, self % latSizeN)

    do i=1, 3
      matBounds(i)   = (ijk(i)-1) * self % latPitch(i) + self % bounds(i) + SURF_TOL
      matBounds(i+3) = ijk(i)     * self % latPitch(i) + self % bounds(i) - SURF_TOL
    end do

  end function getMatBounds

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(materialSource), intent(inout) :: self

    call kill_super(self)

    self % isMG   = .true.
    self % bounds = ZERO
    self % G      = 0

  end subroutine kill

  !!
  !! Generate ijk from flipped localID and shape
  !!
  !! Note that this is NOT the same as get_ijk in latUniverse_class. Lattice is built with first
  !! map input as x_min, y_MAX, z_MAX cell, but localID begins at x_min, y_min, z_min cell. In
  !! this module we want to find ijk from matIdx, which we convert to a flippedLocalID by
  !! offsetting for void regions, which starts counting from the wrong corner. We therefore flip
  !! ijk in the y and z directions in this function compared to instances of this function in other
  !! modules.
  !!
  !! Args:
  !!   flippedlocalID [in] -> Local id of the cell between 1 and product(sizeN),
  !!                          counting from wrong corner
  !!   sizeN [in]          -> Number of cells in each cardinal direction x, y & z
  !!
  !! Result:
  !!   Array ijk which has integer position in each cardinal direction
  !!
  pure function get_ijk(flippedLocalID, sizeN) result(ijk)
    integer(shortInt), intent(in)               :: flippedLocalID
    integer(shortInt), dimension(3), intent(in) :: sizeN
    integer(shortInt), dimension(3)             :: ijk
    integer(shortInt)                           :: temp, base

    temp = flippedLocalID - 1
    base = temp / sizeN(1)
    ijk(1) = temp - sizeN(1) * base + 1

    temp = base
    base = temp / sizeN(2)
    ijk(2) = sizeN(2)*(1 + base) - temp

    ijk(3) = sizeN(3) - base

  end function get_ijk

end module materialSource_class
