module bbSurfaceSource_class

  use numPrecision
  use universalVariables
  use genericProcedures,     only : fatalError, numToChar
  use particle_class,        only : particleState, P_NEUTRON, P_PHOTON
  use particleDungeon_class, only : particleDungeon
  use dictionary_class,      only : dictionary
  use configSource_inter,    only : configSource, kill_super => kill
  use geometry_inter,        only : geometry
  use RNG_class,             only : RNG

  implicit none
  private

  !!
  !! Generates a source representing a black body surface
  !!
  !! Private members:
  !!   r            -> bottom corner of source
  !!   dr           -> size of surface, will be 0 in one dimension 
  !!   dir          -> direction of dominant movement: [1,0,0], [-1,0,0], [0,1,0], etc.
  !!   particleType -> source particle type (photon)
  !!   isMG         -> is the source multi-group? (yes)
  !!
  !! Interface:
  !!   init              -> initialise point source
  !!   append            -> source particles and add to existing dungeon
  !!   sampleType        -> set particle type
  !!   samplePosition    -> set particle position
  !!   sampleEnergyAngle -> sample particle angle
  !!   sampleEnergy      -> set particle energy (isMG = .true., G = 1)
  !!   sampleWeight      -> set particle energy-weight
  !!   kill              -> terminate source
  !!
  !! Sample Dictionary Input:
  !!   source {
  !!       type bbSurfaceSource;
  !!       r (x_min x_max y_min y_max z_min z_max); -> Position bounds of surface
  !!                     -> min and max must be equal in one dimension
  !!       #dir -1;      -> optional, negative will reverse direction in dominant axis
  !!                     -> defaults to positive
  !!       temp 1;       -> temperature of the black body source
  !!       #deltaT 0.05; -> time step size, automatically added to dictionary in IMCPhysicsPackage_class.f90
  !!       N 100;        -> number of particles per time step, only used if append is called with N = 0
  !!      }
  !!
  type, public,extends(configSource) :: bbSurfaceSource
    private
    real(defReal), dimension(3)     :: r            = ZERO
    real(defReal), dimension(3)     :: dr           = ZERO
    integer(shortInt), dimension(3) :: dir          = ZERO
    integer(shortInt)               :: particleType = P_PHOTON
    logical(defBool)                :: isMG         = .true.
    real(defReal)                   :: T            = ZERO
    real(defReal)                   :: deltaT       = ZERO
    integer(shortInt)               :: N            = 0
  contains
    procedure :: init
    procedure :: append
    procedure :: sampleType
    procedure :: samplePosition
    procedure :: sampleEnergy
    procedure :: sampleEnergyAngle
    procedure :: sampleWeight
    procedure :: kill
  end type bbSurfaceSource

contains

  !!
  !! Initialise from dictionary
  !!
  !! See source_inter for details
  !!
  !! Errors:
  !!   - error if an unrecognised particle type is provided
  !!   - error if an axis other than x, y, or z is given
  !!   - error if shape is not square or circle
  !!
  subroutine init(self, dict, geom)
    class(bbSurfaceSource), intent(inout)    :: self
    class(dictionary), intent(in)            :: dict
    class(geometry), pointer, intent(in)     :: geom
    character(30)                            :: type, tempName
    integer(shortInt)                        :: matIdx, uniqueID
    logical(defBool)                         :: isCE, isMG
    real(defReal), dimension(:), allocatable :: temp
    integer(shortInt)                        :: i, dir
    character(100), parameter :: Here = 'init (bbSurfaceSource_class.f90)'

    ! Provide geometry info to source
    self % geom => geom

    ! Provide particle type
    self % particleType = P_PHOTON

    ! Get and check position vector
    call dict % get(temp, 'r')
    if (size(temp) /= 6) call fatalError(Here, 'r should be of size 6')
    do i = 1, 3
      ! Store x_min, y_min, z_min
      self % r(i) = temp(2*i-1)
      ! Store dx, dy, dz
      self % dr(i) = temp(2*i) - temp(2*i-1)
      ! Check for compatible min and max
      if (self % dr(i) < 0) call fatalError(Here, 'Min > Max along direction '//numToChar(i))
    end do
    ! Check that exactly one normal axis is present
    if (count(self % dr == 0) /= 1) call fatalError(Here, 'No clearly defined axis extracted')

    ! Get primary direction
    call dict % getOrDefault(dir, 'dir', 1)
    do i = 1, 3
      if (self % dr(i) == 0) self % dir(i) = sign(1, dir)
    end do

    ! Move by 2*SURF_TOL to ensure sourcing in correct material
    self % r = self % r + 2*SURF_TOL*self % dir

    ! Get remaining information
    call dict % get(self % T, 'temp')
    call dict % get(self % deltaT, 'deltaT') ! Automatically added to dict in IMC physics package
    call dict % getOrDefault(self % N, 'N', 1)

  end subroutine init

  !!
  !! Add particles to given dungeon
  !!
  !! See source_inter for details
  !!
  !! If N is given as 0, then N is instead taken from the input dictionary defining this source
  !! to allow PP to have control over particle numbers
  !!
  subroutine append(self, dungeon, N, rand)
    class(bbSurfaceSource), intent(inout)   :: self
    type(particleDungeon), intent(inout)    :: dungeon
    integer(shortInt), intent(in)           :: N
    class(RNG), intent(inout)               :: rand
    integer(shortInt)                       :: i
    type(RNG)                               :: pRand
    character(100), parameter               :: Here = 'append (bbSurfaceSource_class.f90)'

    ! Set number to generate. Using 0 in function call will use N from input dictionary
    if (N /= 0) self % N = N
    ! TODO change so that this override is only temporary, so that can be called with 0 again later

    ! Generate N particles to populate dungeon
    !$omp parallel
    pRand = rand
    !$omp do private(pRand)
    do i = 1, self % N
      call pRand % stride(i)
      call dungeon % detain(self % sampleParticle(pRand))
    end do
    !$omp end do 
    !$omp end parallel

  end subroutine append

  !!
  !! Provide particle type
  !!
  !! See configSource_inter for details.
  !!
  subroutine sampleType(self, p, rand)
    class(bbSurfaceSource), intent(inout)   :: self
    class(particleState), intent(inout) :: p
    class(RNG), intent(inout)           :: rand

    p % type = self % particleType

  end subroutine sampleType

  !!
  !! Provide particle position
  !!
  !! See configSource_inter for details.
  !!
  subroutine samplePosition(self, p, rand)
    class(bbSurfaceSource), intent(inout) :: self
    class(particleState), intent(inout)   :: p
    class(RNG), intent(inout)             :: rand
    integer(shortInt)                     :: i
    real(defReal), dimension(3)           :: r

    ! Set new x, y and z coords
    do i = 1, 3
      r(i) = (self % dr(i)) * rand % get() + self % r(i)
    end do

    ! Assign to particle
    p % r = r

  end subroutine samplePosition

  !!
  !! Sample angle
  !!
  !! See configSource_inter for details.
  !!
  subroutine sampleEnergyAngle(self, p, rand)
    class(bbSurfaceSource), intent(inout) :: self
    class(particleState), intent(inout)   :: p
    class(RNG), intent(inout)             :: rand
    real(defReal), dimension(3)           :: dir
    real(defReal)                         :: phi, mu
    character(100), parameter :: Here = 'sampleEnergyAngle (bbSurfaceSource_class.f90)'

    ! Sample required phi and mu
    phi = TWO_PI * rand % get()
    mu = sqrt(rand % get())

    ! Choose direction based on dominant direction given in self % dir
    if      (self % dir(1) ==  1) then  ! Positive x
      dir = [ mu, sqrt(1-mu*mu)*cos(phi), sqrt(1-mu*mu)*sin(phi)]

    else if (self % dir(1) == -1) then  ! Negative x
      dir = [-mu, sqrt(1-mu*mu)*cos(phi), sqrt(1-mu*mu)*sin(phi)]

    else if (self % dir(2) ==  1) then  ! Positive y
      dir = [sqrt(1-mu*mu)*sin(phi),  mu, sqrt(1-mu*mu)*cos(phi)]

    else if (self % dir(2) == -1) then  ! Negative y
      dir = [sqrt(1-mu*mu)*sin(phi), -mu, sqrt(1-mu*mu)*cos(phi)]

    else if (self % dir(3) ==  1) then  ! Positive z
      dir = [sqrt(1-mu*mu)*cos(phi), sqrt(1-mu*mu)*sin(phi),  mu]

    else if (self % dir(3) == -1) then  ! Negative z
      dir = [sqrt(1-mu*mu)*cos(phi), sqrt(1-mu*mu)*sin(phi), -mu]

    else
      call fatalError(Here, 'Invalid direction vector')
    end if

    ! Assign to particle
    p % dir = dir

  end subroutine sampleEnergyAngle

  !!
  !! Provide particle energy, currently only a single group
  !!
  !! See configSource_inter for details.
  !!
  subroutine sampleEnergy(self, p, rand)
    class(bbSurfaceSource), intent(inout) :: self
    class(particleState), intent(inout)   :: p
    class(RNG), intent(inout)             :: rand
    real(defReal)                         :: num

    p % isMG = .true.
    p % G    = 1

  end subroutine sampleEnergy

  !!
  !! Provide particle energy-weight
  !!
  !! Sampled as a black body surface, see "Four Decades of Implicit Monte Carlo",
  !!  Allan B Wollaber, p.24-25
  !!
  !! See configSource_inter for details.
  !!
  subroutine sampleWeight(self, p, rand)
    class(bbSurfaceSource), intent(inout) :: self
    class(particleState), intent(inout)   :: p
    class(RNG), intent(inout)             :: rand
    real(defReal)                         :: area, num

    ! Calculate surface area of source
    area = product(self % dr, self % dr /= ZERO)

    ! Calculate energy weight per particle
    num = radiationConstant * lightSpeed * self % deltaT * self % T**4 * area
    p % wgt = num / (4 * self % N)

  end subroutine sampleWeight

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(bbSurfaceSource), intent(inout) :: self

    ! Kill superclass
    call kill_super(self)

    ! Kill local components
    self % r   = ZERO
    self % dr  = ZERO
    self % dir = ZERO
    self % particleType = P_PHOTON
    self % isMG         = .true.
    self % T      = ZERO
    self % deltaT = ZERO
    self % N      = ZERO

  end subroutine kill

end module bbSurfaceSource_class
