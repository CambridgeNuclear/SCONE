module blackBodySource_class

  use numPrecision
  use universalVariables
  use genericProcedures,     only : fatalError, numToChar
  use particle_class,        only : particleState, P_NEUTRON, P_PHOTON
  use particleDungeon_class, only : particleDungeon
  use dictionary_class,      only : dictionary
  use configSource_inter,    only : configSource, kill_super => kill
  use geometry_inter,        only : geometry
  use RNG_class,             only : RNG
  use simulationTime_class

  use energyGrid_class,         only : energyGrid
  use energyGridRegistry_mod,   only : get_energyGrid

  implicit none
  private

  ! Options for source distribution
  integer(shortInt), parameter, public :: SURFACE = 1
  integer(shortInt), parameter, public :: UNIFORM = 2
  integer(shortInt), parameter, public :: OLSON1D = 3

  !!
  !! Generates a source representing a black body
  !! Standard is a surface distribution but can be configured for custom distributions
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
  !!   sampleEnergy      -> set particle energy
  !!   sampleWeight      -> set particle energy-weight
  !!   sampleTime        -> set particle time
  !!   kill              -> terminate source
  !!
  !! Sample Dictionary Input:
  !!   source {
  !!       type blackBodySource;
  !!       distribution surface;
  !!       surface -x;
  !!       temp 1;       -> temperature of the black body source
  !!      }
  !!
  !!  Current source distributions:
  !!    surface -> black body surface source placed on the surface indicated
  !!    uniform -> uniform in space, isotropic
  !!    olson1D -> 1D multifrequency benchmark source from "Stretched and Filtered Multigroup Pn
  !!               Transport for Improved Positivity and Accuracy", Olson, Gordon Lee, 2020
  !!               See materialEquations.f90 for sigma and cv equations for this benchmark.
  !!
  type, public,extends(configSource) :: blackBodySource
    private

    ! Settings defining sourced position and direction for surface sources
    real(defReal), dimension(3)       :: r            = ZERO     ! Corner position of source
    real(defReal), dimension(3)       :: dr           = ZERO     ! Spatial extent from corner
    integer(shortInt), dimension(3,3) :: rotation     = ZERO     ! Direction rotation matrix
    ! Other settings
    integer(shortInt)                 :: distribution = SURFACE
    integer(shortInt)                 :: particleType = P_PHOTON
    logical(defBool)                  :: isMG         = .true.
    real(defReal)                     :: T            = ZERO     ! Source temperature
    real(defReal)                     :: particleWeight = ZERO   ! Weight of each particle (= sourceWeight/N)
    integer(shortInt)                 :: timeStepMax  = 0        ! Time step to switch source off
    ! Probability of emission from each energy group for multi-frequency case
    real(defReal), dimension(:), allocatable :: cumEnergyProbs

  contains
    ! Initialisation procedures
    procedure :: init
    procedure :: initSurface
    procedure :: initCustom
    ! Sampling procedures
    procedure :: append
    procedure :: sampleType
    procedure :: samplePosition
    procedure :: sampleEnergy
    procedure :: sampleEnergyAngle
    procedure :: sampleWeight
    procedure :: sampleTime
    procedure :: kill

  end type blackBodySource

contains

!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Particle sampling procedures
!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Add particles to given dungeon
  !!
  !! See source_inter for details
  !!
  !! If N is given as 0, then N is instead taken from the input dictionary defining this source
  !! to allow PP to have control over particle numbers
  !!
  subroutine append(self, dungeon, N, rand)
    class(blackBodySource), intent(inout)   :: self
    type(particleDungeon), intent(inout)    :: dungeon
    integer(shortInt), intent(in)           :: N
    class(RNG), intent(inout)               :: rand
    integer(shortInt)                       :: i
    type(RNG)                               :: pRand
    character(100), parameter               :: Here = 'append (blackBodySource_class.f90)'

    ! Store N for calculation of each particle weight
    self % particleWeight = self % sourceWeight / N

    ! Generate N particles to populate dungeon
    !$omp parallel
    pRand = rand
    !$omp do private(pRand)
    do i = 1, N
      call pRand % stride(i)
      call dungeon % detain(self % sampleParticle(pRand))
    end do
    !$omp end do 
    !$omp end parallel

    ! Turn off for next time step if needed
    if (thisStep() == self % timeStepMax) self % sourceWeight = ZERO

  end subroutine append

  !!
  !! Provide particle type
  !!
  !! See configSource_inter for details.
  !!
  subroutine sampleType(self, p, rand)
    class(blackBodySource), intent(inout) :: self
    class(particleState), intent(inout)   :: p
    class(RNG), intent(inout)             :: rand

    p % type = self % particleType

  end subroutine sampleType

  !!
  !! Provide particle position
  !!
  !! See configSource_inter for details.
  !!
  subroutine samplePosition(self, p, rand)
    class(blackBodySource), intent(inout) :: self
    class(particleState), intent(inout)   :: p
    class(RNG), intent(inout)             :: rand
    integer(shortInt)                     :: i
    real(defReal), dimension(3)           :: r
    real(defReal)                         :: x
    character(100), parameter :: Here = 'samplePosition (blackBodySource_class.f90)'

    select case(self % distribution)

      case(SURFACE, UNIFORM)
        ! Set new x, y and z coords
        do i = 1, 3
          r(i) = self % r(i) + rand % get()*self % dr(i)
        end do

      case(OLSON1D)
        ! Q(x) proportional to exp(-693x**3) (integral from 0 to 4.8 = 0.100909)
        rejection:do
          x = rand % get() * 4.8_defReal
          if (rand % get() < exp(-693*x**3)/0.100909) exit
        end do rejection
        r(1) = x - 2.4_defReal
        r(2) = rand % get() - 0.5_defReal
        r(3) = rand % get() - 0.5_defReal

      case default
        call fatalError(Here, 'Unrecognised source distribution')

    end select

    ! Assign to particle
    p % r = r

  end subroutine samplePosition

  !!
  !! Sample angle
  !!
  !! See configSource_inter for details.
  !!
  subroutine sampleEnergyAngle(self, p, rand)
    class(blackBodySource), intent(inout) :: self
    class(particleState), intent(inout)   :: p
    class(RNG), intent(inout)             :: rand
    real(defReal), dimension(3)           :: dir
    real(defReal)                         :: phi, mu
    character(100), parameter :: Here = 'sampleEnergyAngle (blackBodySource_class.f90)'

    ! Sample direction for isotropic source
    if (all(self % rotation == ZERO)) then
      ! Sample uniformly within unit sphere
      mu = 2 * rand % get() - 1
      phi = rand % get() * 2*pi
      dir(1) = mu
      dir(2) = sqrt(1-mu**2) * cos(phi)
      dir(3) = sqrt(1-mu**2) * sin(phi)
      ! Assign to particle
      p % dir = dir
      return
    end if

    ! If not isotropic (e.g. surface distribution), sample first with a primary direction of +x
    phi = TWO_PI * rand % get()
    mu = sqrt(rand % get())
    dir = [mu, sqrt(1-mu*mu)*cos(phi), sqrt(1-mu*mu)*sin(phi)]

    ! Rotate to direction of surface normal
    p % dir = matmul(self % rotation, dir)

  end subroutine sampleEnergyAngle

  !!
  !! Provide particle energy
  !!
  !! See configSource_inter for details.
  !!
  subroutine sampleEnergy(self, p, rand)
    class(blackBodySource), intent(inout) :: self
    class(particleState), intent(inout)   :: p
    class(RNG), intent(inout)             :: rand
    real(defReal)                         :: random
    integer(shortInt)                     :: i
    character(100), parameter :: Here = 'sampleEnergy (blackBodySource_class.f90)'

    p % isMG = .true.

    ! Grey case
    if (.not.allocated(self % cumEnergyProbs)) then
      p % G = 1
      return
    end if

    ! Sample energy group
    random = rand % get()
    do i=1, size(self % cumEnergyProbs)
      if (random <= self % cumEnergyProbs(i)) then
        p % G = i
        return
      end if
    end do

    call fatalError(Here, 'Somehow failed to sample particle energy group')

  end subroutine sampleEnergy

  !!
  !! Provide particle energy-weight
  !!
  !! See configSource_inter for details.
  !!
  subroutine sampleWeight(self, p, rand)
    class(blackBodySource), intent(inout) :: self
    class(particleState), intent(inout)   :: p
    class(RNG), intent(inout)             :: rand

    p % wgt = self % particleWeight

  end subroutine sampleWeight

  !!
  !! Sample time uniformly within time step
  !!
  subroutine sampleTime(self, p, rand)
    class(blackBodySource), intent(inout) :: self
    class(particleState), intent(inout)   :: p
    class(RNG), intent(inout)             :: rand

    ! Sample time uniformly within time step
    p % time = time % stepStart + timeStep() * rand % get()

  end subroutine sampleTime


!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Source initialisation procedures
!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Initialise from dictionary
  !!
  !! See source_inter for details
  !!
  subroutine init(self, dict, geom)
    class(blackBodySource), intent(inout)    :: self
    class(dictionary), intent(in)            :: dict
    class(geometry), pointer, intent(in)     :: geom
    integer(shortInt)                        :: i, j, nGroups
    real(defReal)                            :: nu, eStep
    type(energyGrid)                         :: eGrid
    logical(defBool)                         :: err
    character(nameLen)                       :: gridName, distribution
    character(100), parameter :: Here = 'init (blackBodySource_class.f90)'

    ! Provide geometry info to source
    self % geom => geom

    ! Provide particle type
    self % particleType = P_PHOTON

    ! Get source temperature
    call dict % get(self % T, 'temp')

    ! Initialise specifics of source (e.g. position samlping bounds)
    call dict % getOrDefault(distribution, 'distribution', 'surface')
    select case(distribution)
      case('surface')
        call self % initSurface(dict)
      case default
        call self % initCustom(dict)
    end select

    ! Get time step to turn off source
    call dict % getOrDefault(self % timeStepMax, 'untilStep', 0)

    ! Exit in grey case
    gridName = 'mg'
    call get_energyGrid(eGrid, gridName, err)
    if (err .eqv. .true.) return

    ! Calculate emission probability in each energy group
    nGroups = eGrid % getSize()
    allocate(self % cumEnergyProbs(nGroups))

    do i=1, nGroups
      eStep = (eGrid % bin(i) - eGrid % bin(i+1)) / 1000
      nu = eGrid % bin(i) - eStep/2
      ! Add previous group probability to get cumulative distribution
      if (i > 1) self % cumEnergyProbs(i) = self % cumEnergyProbs(i-1)
      ! Step through energies
      do j=1, 1000
        self % cumEnergyProbs(i) = self % cumEnergyProbs(i) + eStep*normPlanck(nu,self % T)
        nu = nu - eStep
      end do
    end do

    ! Normalise to account for exclusion of tail ends of distribution
    ! TODO: should possibly add these tails into outer energy groups rather than normalising over all groups?
    self % cumEnergyProbs = self % cumEnergyProbs / self % cumEnergyProbs(nGroups)

  end subroutine init


  !!
  !! Initialise source for standard black body surface by placing source as one side of
  !! bounding box of geometry
  !!
  !! Input dict should contain 'surface', corresponding to which side of the box is the source
  !! e.g. surface -x;    => source placed on negative x side of bounding box
  !!      surface +z;    => source placed on positive z side of bounding box
  !!      etc.
  !!
  !! Automatically sets bounds for position sampling, rotation matrix for direction sampling, and
  !! total weight emitted from the source during each time step
  !!
  subroutine initSurface(self, dict)
    class(blackBodySource), intent(inout)   :: self
    class(dictionary), intent(in)           :: dict
    real(defReal), dimension(6)             :: bounds
    character(nameLen)                      :: whichSurface
    integer(shortInt)                       :: i
    integer(shortInt), dimension(9)         :: rotation
    real(defReal)                           :: area
    character(100), parameter :: Here = 'initSurface (blackBodySource_class.f90)'

    self % distribution = SURFACE

    bounds = self % geom % bounds()

    ! Bottom left corner
    self % r = bounds(1:3)
    ! Dimensions of bounding box
    do i = 1, 3
      self % dr(i) = bounds(i+3) - bounds(i)
    end do

    ! Get position bounds
    call dict % get(whichSurface, 'surface')

    select case(whichSurface)
      case('-x')
        ! Set sampling position to be at constant x value
        self % dr(1) = ZERO
        ! Nudge to ensure sourcing in correct material
        self % r(1) = bounds(1) + TWO*SURF_TOL
        ! Set rotation matrix for direction sampling
        rotation = [[1,0,0],[0,1,0],[0,0,1]]
      case('-y')
        self % dr(2) = ZERO
        self % r(2) = bounds(2) + TWO*SURF_TOL
        rotation = [[0,1,0],[1,0,0],[0,0,1]]
      case('-z')
        self % dr(3) = ZERO
        self % r(3) = bounds(3) + TWO*SURF_TOL
        rotation = [[0,0,1],[0,1,0],[1,0,0]]
      case('+x')
        self % r(1)  = bounds(4) - TWO*SURF_TOL
        self % dr(1) = ZERO
        rotation = [[-1,0,0],[0,1,0],[0,0,1]]
      case('+y')
        self % r(2)  = bounds(5) - TWO*SURF_TOL
        self % dr(2) = ZERO
        rotation = [[0,-1,0],[1,0,0],[0,0,1]]
      case('+z')
        self % r(3)  = bounds(6) - TWO*SURF_TOL
        self % dr(3) = ZERO
        rotation = [[0,0,-1],[0,1,0],[1,0,0]]
      case default
        call fatalError(Here, 'Unrecognised surface, must be +/-x,y or z')
    end select

    ! Note that reshape fills columns first so the above rotations are [[col1][col2][col3]]
    self % rotation = reshape(rotation,[3,3])

    ! Calculate surface area of source
    area = product(self % dr, self % dr /= ZERO)

    ! Calculate total source energy
    self % sourceWeight = radiationConstant * lightSpeed * timeStep() * self % T**4 * area / 4

  end subroutine initSurface

  !!
  !! Initialise black body source for more unique cases, e.g. Olson 1D MG benchmark (2020)
  !!
  !! Rotation matrix of ZERO corresponds to isotropic direction sampling. Set source weight or
  !! other settings as required for specific case, and add additional cases to individual samping
  !! procedures if needed.
  !!
  subroutine initCustom(self, dict)
    class(blackBodySource), intent(inout) :: self
    class(dictionary), intent(in)         :: dict
    character(nameLen)                    :: name
    integer(shortInt)                     :: i
    real(defReal), dimension(6)           :: bounds
    character(100), parameter :: Here = 'initCustom (blackBodySource_class.f90)'

    call dict % get(name, 'distribution')

    select case(name)

      case('uniform')
        self % distribution = UNIFORM
        bounds = self % geom % bounds()
        ! Bottom left corner
        self % r = bounds(1:3)
        ! Dimensions of bounding box
        do i = 1, 3
          self % dr(i) = bounds(i+3) - bounds(i)
        end do
        ! Isotropic direction sampling
        self % rotation = ZERO
        call dict % get(self % sourceWeight, 'sourceWeight')

      case('olson1D')
        ! See self % samplePosition for position sampling specifics for Olson1D
        self % distribution = OLSON1D
        ! Isotropic directional sampling
        self % rotation = ZERO
        ! Set source weight
        self % sourceWeight = radiationConstant * lightSpeed * timeStep() * self % T**4 * 0.100909

      case default
        call fatalError(Here, 'Unrecognised distribution for custom black body source')

    end select

  end subroutine initCustom

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(blackBodySource), intent(inout) :: self

    ! Kill superclass
    call kill_super(self)

    ! Kill local components
    self % r              = ZERO
    self % dr             = ZERO
    self % distribution   = SURFACE
    self % particleType   = P_PHOTON
    self % isMG           = .true.
    self % T              = ZERO
    self % particleWeight = ZERO

  end subroutine kill


!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Misc. procedures
!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Return the frequency-normalised Planck spectrum evaluated for given frequency
  !! nu and temperature T
  !!
  pure function normPlanck(E, T) result(b)
    real(defReal), intent(in) :: E, T
    real(defReal)             :: b

    b = 15*E**3 / ((pi*T)**4 * (exp(E/T)-1))

  end function normPlanck


end module blackBodySource_class
