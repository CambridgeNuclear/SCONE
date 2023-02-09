module IMCSource_class

  use numPrecision
  use endfConstants
  use universalVariables
  use genericProcedures,       only : fatalError, rotateVector
  use dictionary_class,        only : dictionary
  use RNG_class,               only : RNG

  use particle_class,          only : particle, particleState, P_PHOTON
  use particleDungeon_class,   only : particleDungeon
  use source_inter,            only : source, kill_super => kill

  use geometry_inter,          only : geometry
  use IMCMaterial_inter,       only : IMCMaterial, IMCMaterial_CptrCast
  use nuclearDataReg_mod,      only : ndReg_getIMCMG => getIMCMG
  use nuclearDatabase_inter,   only : nuclearDatabase
  use mgIMCDatabase_inter,     only : mgIMCDatabase

  implicit none
  private

  !!
  !! IMC Source for uniform generation of photons within a material
  !!
  !! Angular distribution is isotropic.
  !!
  !! Private members:
  !!   isMG    -> is the source multi-group? (default = .true.)
  !!   bottom  -> Bottom corner (x_min, y_min, z_min)
  !!   top     -> Top corner (x_max, y_max, z_max)
  !!   G       -> Group (default = 1)
  !!   N       -> number of particles being generated, used to normalise weight in sampleParticle
  !!   matIdx  -> index of material to be sampled from
  !!
  !! Interface:
  !!   source_inter Interface
  !!
  !! SAMPLE INPUT:
  !!   imcSource { type IMCSource; }
  !!
  type, public,extends(source) :: imcSource
    private
    logical(defBool)                             :: isMG   = .true.
    real(defReal), dimension(3)                  :: bottom = ZERO
    real(defReal), dimension(3)                  :: top    = ZERO
    real(defReal), dimension(3)                  :: latPitch = ZERO
    integer(shortInt), dimension(:), allocatable :: latSizeN
    integer(shortInt)                            :: G      = 0
    integer(shortInt)                            :: N
    integer(shortInt)                            :: matIdx
    real(defReal), dimension(6)                  :: matBounds = ZERO
  contains
    procedure :: init
    procedure :: append
    procedure :: sampleParticle
    procedure :: samplePosRej
    procedure :: samplePosLat
    procedure :: kill
  end type imcSource

contains

  !!
  !! Initialise IMC Source
  !!
  !! See source_inter for details
  !!
  subroutine init(self, dict, geom)
    class(imcSource), intent(inout)          :: self
    class(dictionary), intent(in)            :: dict
    class(geometry), pointer, intent(in)     :: geom
    real(defReal), dimension(6)              :: bounds
    character(100), parameter :: Here = 'init (imcSource_class.f90)'

    ! Provide geometry info to source
    self % geom => geom

    call dict % getOrDefault(self % G, 'G', 1)

    ! Set bounding region
    bounds = self % geom % bounds()
    self % bottom = bounds(1:3)
    self % top    = bounds(4:6)

    ! Store lattice dimensions for use in position sampling if using a large lattice
    ! sizeN automatically added to dict in IMCPhysicsPackage if needed
    if (dict % isPresent('sizeN')) then
      call dict % get(self % latSizeN, 'sizeN')
      self % latPitch = (self % top - self % bottom) / self % latSizeN
    end if

  end subroutine init

  !!
  !! Generate N particles from material matIdx to add to a particleDungeon without overriding
  !! particles already present.
  !!
  !! Args:
  !!   dungeon [inout] -> particle dungeon to be added to
  !!   n [in]          -> number of particles to place in dungeon
  !!   rand [inout]    -> particle RNG object
  !!   matIdx [in]     -> index of material to sample from
  !!
  !! Result:
  !!   A dungeon populated with N particles sampled from the source, plus particles
  !!   already present in dungeon
  !!
  subroutine append(self, dungeon, N, rand, matIdx)
    class(imcSource), intent(inout)         :: self
    type(particleDungeon), intent(inout)    :: dungeon
    integer(shortInt), intent(in)           :: N
    class(RNG), intent(inout)               :: rand
    integer(shortInt), intent(in), optional :: matIdx
    integer(shortInt)                       :: i
    integer(shortInt), dimension(3)         :: ijk
    type(RNG)                               :: pRand
    character(100), parameter               :: Here = "append (IMCSource_class.f90)"

    ! Assert that optional argument matIdx is in fact present
    if (.not. present(matIdx)) call fatalError(Here, 'matIdx must be provided for IMC source')

    ! Store inputs for use by sampleParticle subroutine
    self % N      = N
    self % matIdx = matIdx

    ! For a large number of materials (large lattice using discretiseGeom_class) rejection
    ! sampling is too slow, so calculate bounding box of material
    if (self % latPitch(1) /= 0) then
      ijk = get_ijk(matIdx, self % latSizeN)
      do i=1, 3
        self % matBounds(i)   = (ijk(i)-1) * self % latPitch(i) + self % bottom(i)
        self % matBounds(i+3) = ijk(i)     * self % latPitch(i) + self % bottom(i)
      end do
    end if

    ! Add N particles to dungeon
    !$omp parallel
    pRand = rand
    !$omp do private(pRand)
    do i=1, N
      call pRand % stride(i)
      call dungeon % detain(self % sampleParticle(pRand))
    end do
    !$omp end do
    !$omp end parallel
   
  end subroutine append

  !!
  !! Sample particle's phase space co-ordinates
  !!
  !! See source_inter for details
  !!
  function sampleParticle(self, rand) result(p)
    class(imcSource), intent(inout)      :: self
    class(RNG), intent(inout)            :: rand
    type(particleState)                  :: p
    class(nuclearDatabase), pointer      :: nucData
    class(IMCMaterial), pointer          :: mat
    real(defReal), dimension(3)          :: r, dir
    real(defReal)                        :: mu, phi
    integer(shortInt)                    :: matIdx, uniqueID
    character(100), parameter :: Here = 'sampleParticle (imcSource_class.f90)'

    ! Get pointer to appropriate nuclear database
    nucData => ndReg_getIMCMG()
    if(.not.associated(nucData)) call fatalError(Here, 'Failed to retrieve Nuclear Database')

    ! Choose position sampling method
    if (self % latPitch(1) == ZERO) then
      call self % samplePosRej(r, matIdx, uniqueID, rand)
    else
      call self % samplePosLat(r, matIdx, uniqueID, rand)
    end if

    ! Point to material
    mat => IMCMaterial_CptrCast(nucData % getMaterial(matIdx))
    if (.not.associated(mat)) call fatalError(Here, "Nuclear data did not return IMC material.")

    ! Sample direction - chosen uniformly inside unit sphere
    mu = 2 * rand % get() - 1
    phi = rand % get() * 2*pi
    dir(1) = mu
    dir(2) = sqrt(1-mu**2) * cos(phi)
    dir(3) = sqrt(1-mu**2) * sin(phi)

    ! Assign basic phase-space coordinates
    p % matIdx   = matIdx
    p % uniqueID = uniqueID
    p % time     = ZERO
    p % type     = P_PHOTON
    p % r        = r
    p % dir      = dir
    p % G        = self % G
    p % isMG     = .true.

    ! Set weight
    p % wgt = mat % getEmittedRad() / self % N

  end function sampleParticle


  !!
  !! Position is sampled by taking a random point from within geometry bounding box
  !! If in correct material, position is accepted
  !!
  subroutine samplePosRej(self, r, matIdx, uniqueID, rand)
    class(imcSource), intent(inout)          :: self
    real(defReal), dimension(3), intent(out) :: r
    integer(shortInt), intent(out)           :: matIdx
    integer(shortInt), intent(out)           :: uniqueID
    class(RNG), intent(inout)                :: rand
    integer(shortInt)                        :: i
    real(defReal), dimension(3)              :: rand3
    character(100), parameter :: Here = 'samplePosRej (IMCSource_class.f90)'

    i = 0

    rejectionLoop : do

      ! Protect against infinite loop
      i = i+1
      if (i > 10000) then
        call fatalError(Here, '10,000 failed samples in rejection sampling loop')
      end if

      ! Sample Position
      rand3(1) = rand % get()
      rand3(2) = rand % get()
      rand3(3) = rand % get()
      r = (self % top - self % bottom) * rand3 + self % bottom

      ! Find material under position
      call self % geom % whatIsAt(matIdx, uniqueID, r)

      ! Exit if in desired material
      if (matIdx == self % matIdx) exit rejectionLoop

    end do rejectionLoop

  end subroutine samplePosRej

  !!
  !! Sample position without using a rejection sampling method, by calculating the material bounds.
  !!
  !! Requires geometry to be a uniform lattice, so currently only called when discretiseGeom_class
  !! is used to create inputs.
  !!
  subroutine samplePosLat(self, r, matIdx, uniqueID, rand)
    class(imcSource), intent(inout)          :: self
    real(defReal), dimension(3), intent(out) :: r
    integer(shortInt), intent(out)           :: matIdx
    integer(shortInt), intent(out)           :: uniqueID
    class(RNG), intent(inout)                :: rand
    integer(shortInt)                        :: i
    character(100), parameter :: Here = 'samplePosLat (IMCSource_class.f90)'

    do i=1, 3
      r(i) = self % matBounds(i) + rand % get() * (self % matBounds(i+3) - self % matBounds(i) - SURF_TOL) + SURF_TOL
    end do

    call self % geom % whatIsAt(matIdx, uniqueID, r)

    if (matIdx /= self % matIdx) call fatalError(Here, 'Incorrect material')

  end subroutine samplePosLat

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(imcSource), intent(inout) :: self

    call kill_super(self)

    self % isMG   = .true.
    self % bottom = ZERO
    self % top    = ZERO
    self % G      = 0

  end subroutine kill


  !!
  !! Generate ijk from localID and shape
  !!
  !! Args:
  !!   localID [in] -> Local id of the cell between 1 and product(sizeN)
  !!   sizeN [in]   -> Number of cells in each cardinal direction x, y & z
  !!
  !! Result:
  !!   Array ijk which has integer position in each cardinal direction
  !!
  pure function get_ijk(localID, sizeN) result(ijk)
    integer(shortInt), intent(in)               :: localID
    integer(shortInt), dimension(3), intent(in) :: sizeN
    integer(shortInt), dimension(3)             :: ijk
    integer(shortInt)                           :: temp, base

    temp = localID - 1

    base = temp / sizeN(1)
    ijk(1) = temp - sizeN(1) * base + 1

    temp = base
    base = temp / sizeN(2)
    ijk(2) = temp - sizeN(2) * base + 1

    ijk(3) = base + 1

  end function get_ijk


end module IMCSource_class
