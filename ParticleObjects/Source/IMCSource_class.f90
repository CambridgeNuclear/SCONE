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
    integer(shortInt)                            :: G      = 0
    integer(shortInt)                            :: N
    integer(shortInt)                            :: matIdx
  contains
    procedure :: init
    procedure :: append
    procedure :: sampleParticle
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
    integer(shortInt)                        :: i
    character(100), parameter :: Here = 'init (imcSource_class.f90)'

    ! Provide geometry info to source
    self % geom => geom

    call dict % getOrDefault(self % G, 'G', 1)

    ! Set bounding region
    bounds = self % geom % bounds()
    self % bottom = bounds(1:3)
    self % top    = bounds(4:6)

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
    type(particle)                          :: p
    integer(shortInt)                       :: i
    real(defReal)                           :: normFactor
    type(RNG)                               :: pRand
    character(100), parameter               :: Here = "append (IMCSource_class.f90)"

    ! Assert that optional argument matIdx is in fact present
    if (.not. present(matIdx)) call fatalError(Here, 'matIdx must be provided for IMC source')

    ! Store inputs for use by sampleParticle subroutine
    self % N      = N
    self % matIdx = matIdx

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
    real(defReal), dimension(3)          :: r, rand3, dir
    real(defReal)                        :: mu, phi
    integer(shortInt)                    :: i, matIdx, uniqueID
    character(100), parameter :: Here = 'sampleParticle (imcSource_class.f90)'

    ! Get pointer to appropriate nuclear database
    nucData => ndReg_getIMCMG()
    if(.not.associated(nucData)) call fatalError(Here, 'Failed to retrieve Nuclear Database')

    ! Position is sampled by taking a random point from within geometry bounding box
    ! If in correct material, position is accepted
    i = 0

    rejection : do

      ! Protect against infinite loop
      i = i+1
      if (i > 100000) then
        call fatalError(Here, '100,000 failed samples in rejection sampling loop')
      end if

      ! Sample Position
      rand3(1) = rand % get()
      rand3(2) = rand % get()
      rand3(3) = rand % get()
      r = (self % top - self % bottom) * rand3 + self % bottom

      ! Find material under position
      call self % geom % whatIsAt(matIdx, uniqueID, r)

      ! Reject if not in desired material
      if (matIdx /= self % matIdx) cycle rejection

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

      ! Exit the loop
      exit rejection

    end do rejection

  end function sampleParticle

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

end module IMCSource_class
