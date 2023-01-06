module bbSurfaceSource_class

  use numPrecision
  use universalVariables
  use genericProcedures,     only : fatalError
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
  !! Put together quite quickly so very specific in use and not perfect
  !!   - Currently only allows a circle or square aligned on x y or z axis, with
  !!     a certain radius or side length
  !!   - May still contain unnecessary lines of code copied from pointSource_class.f90
  !!
  !! Private members:
  !!   r            -> source position
  !!   dir          -> optional source direction
  !!   particleType -> source particle type
  !!   isMG         -> is the source multi-group?
  !!   isIsotropic  -> is the source isotropic?
  !!
  !! Interface:
  !!   init              -> initialise point source
  !!   sampleType        -> set particle type
  !!   samplePosition    -> set particle position
  !!   sampleEnergy      -> set particle energy
  !!   sampleEnergyAngle -> sample particle angle
  !!   kill              -> terminate source
  !!
  !! Sample Dictionary Input:
  !!   source {
  !!       type bbSurfaceSource;
  !!       shape circle    ! circle or square;
  !!       size 5;         ! radius(circle) or side length(square)
  !!       axis x;         ! axis normal to planar shape
  !!       pos 0;          ! distance along axis to place plane
  !!       T 1;            ! temperature of source boundary
  !!       particle photon;
  !!       # dir 1; #      ! Positive or negative to indicate direction along axis
  !!                         If 0 then emit in both directions
  !!       # N 100; #      ! Number of particles, only used if call to append subroutine uses N=0
  !!      }
  !!
  type, public,extends(configSource) :: bbSurfaceSource
    private
    real(defReal),dimension(3)  :: r            = ZERO
    real(defReal)               :: dir          = ZERO
    real(defReal)               :: surfSize     = ZERO
    real(defReal)               :: area         = ZERO
    integer(shortInt)           :: particleType = P_PHOTON
    logical(defBool)            :: isMG         = .true.
    logical(defBool)            :: isIsotropic  = .false.
    integer(shortInt)           :: planeShape   = 0  ! 0 => square, 1 => circle
    integer(shortInt)           :: axis         = 1  ! 1 => x, 2 => y, 3 => z
    real(defReal)               :: T            = ZERO
    real(defReal)               :: deltaT       = ZERO
    integer(shortInt)           :: N            = 1
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
    class(bbSurfaceSource), intent(inout) :: self
    class(dictionary), intent(in)         :: dict
    class(geometry), pointer, intent(in)  :: geom
    character(30)                         :: type, tempName
    integer(shortInt)                     :: matIdx, uniqueID
    logical(defBool)                      :: isCE, isMG
    real(defReal)                         :: temp
    character(100), parameter :: Here = 'init (bbSurfaceSource_class.f90)'

    ! Provide geometry info to source
    self % geom => geom

    ! Identify which particle is used in the source
    ! Presently limited to neutron and photon
    call dict % getOrDefault(type, 'particle' ,'photon')
    select case(type)
      case('neutron')
        self % particleType = P_NEUTRON

      case('photon')
        self % particleType = P_PHOTON

      case default
        call fatalError(Here, 'Unrecognised particle type')

    end select

    ! Get position of surface along axis
    call dict % get(temp, 'pos')

    ! Get axis and assign axis position
    call dict % getOrDefault(tempName, 'axis', 'x')
    select case(tempName)
      case('x')
        self % r(1) = temp
        self % axis = 1
      case('y')
        self % r(2) = temp
        self % axis = 2
      case('z')
        self % r(3) = temp
        self % axis = 3
      case default
        call fatalError(Here, 'Unrecognised axis, may only be x, y or z')
    end select

    ! Get size of boundary surface
    call dict % get(self % surfSize, 'size')

    ! Get shape and area of boundary surface
    call dict % get(tempName, 'shape')
    if (tempName == 'square') then
      self % planeShape = 0
      self % area = self % surfSize**2
    else if (tempName == 'circle') then
      self % planeShape = 1
      self % area = pi * self % surfSize**2
    else
      call fatalError(Here, 'Shape must be "square" or "circle"')
    end if

    ! Determine if dir is positive or negative along given axis
    ! If equal to 0, emit from both sides
    self % isIsotropic = .not. dict % isPresent('dir')
    if (.not. self % isIsotropic) then

      call dict % get(temp, 'dir')

      if (temp == 0) then
        self % dir = 0
      else
        ! Set equal to +1 or -1
        self % dir = temp/abs(temp)
      end if

    end if

    call dict % get(self % T, 'T')
    call dict % get(self % deltaT, 'deltaT')
    call dict % getOrDefault(self % N, 'N', 1)

  end subroutine init

  !!
  !! Add particles to given dungeon
  !!
  !! See source_inter for details
  !!
  !! If N is given as 0, then N is instead taken from the input dictionary defining this source
  !!
  subroutine append(self, dungeon, N, rand, matIdx)
    class(bbSurfaceSource), intent(inout)   :: self
    type(particleDungeon), intent(inout)    :: dungeon
    integer(shortInt), intent(in)           :: N
    class(RNG), intent(inout)               :: rand
    integer(shortInt), intent(in), optional :: matIdx
    integer(shortInt)                       :: i
    type(RNG)                               :: pRand
    character(100), parameter               :: Here = 'append (bbSurfaceSource_class.f90)'

    ! Set number to generate. Using 0 in function call will use N from input dictionary
    if (N /= 0) self % N = N


! TODO Parallel for some reason isn't working here, even though changes are the same as IMCSource ???

    ! Generate N particles to populate dungeon
!    !$omp parallel
!    pRand = rand
!    !$omp do private(pRand)
!    do i = 1, self % N
!      call pRand % stride(i)
!      call dungeon % detain(self % sampleParticle(pRand))
!    end do
!    !$omp end do 
!    !$omp end parallel


    do i = 1, self % N
      call dungeon % detain(self % sampleParticle(rand))
    end do

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
    real(defReal), dimension(3)           :: prevPos
    real(defReal)                         :: r1, r2, rad, theta

    if ( self % planeShape == 0 ) then   ! Square

      prevPos = self % r

      ! Set new x, y and z coords
      self % r(1) = (rand % get()-0.5) * self % surfSize
      self % r(2) = (rand % get()-0.5) * self % surfSize
      self % r(3) = (rand % get()-0.5) * self % surfSize
      ! Leave position along normal axis unchanged
      self % r(self % axis) = prevPos(self % axis) 

    else   ! Circle
      rad = rand % get() * self % surfSize
      theta = rand % get() * 2 * pi

      r1 = rad * cos(theta)
      r2 = rad * sin(theta)

      if(self % axis == 1) then  ! Set y and z
        self % r(2) = r1
        self % r(3) = r2
      else if(self % axis == 2) then  ! Set x and z
        self % r(1) = r1
        self % r(3) = r2
      else  ! Set x and y
        self % r(1) = r1
        self % r(2) = r2
      end if

    end if

    p % r = self % r

  end subroutine samplePosition

  !!
  !! Provide angle or sample if isotropic
  !!
  !! See configSource_inter for details.
  !!
  !! Only isotropic/fixed direction. Does not sample energy.
  !!
  subroutine sampleEnergyAngle(self, p, rand)
    class(bbSurfaceSource), intent(inout) :: self
    class(particleState), intent(inout)   :: p
    class(RNG), intent(inout)             :: rand
    real(defReal)                         :: phi, mu

    phi = TWO_PI * rand % get()
    mu = sqrt(rand % get())

    p % dir = [mu, sqrt(1-mu**2)*cos(phi), sqrt(1-mu**2)*sin(phi)]

    ! If dir not equal to zero, adjust so that particles are travelling in correct direction
    if (self % dir /= 0) then
      p % dir(self % axis) = abs(p % dir(self % axis)) * self % dir 
    end if


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
    real(defReal)                         :: num

    num = radiationConstant * lightSpeed * self % deltaT * self % T**4 * self % area
    p % wgt = num / (4 * self % N)

    ! If dir = 0 then emit in both directions => double total energy
    if (self % dir == 0) p % wgt = 2*p % wgt

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
    self % dir = ZERO
    self % particleType = P_PHOTON
    self % isMG         = .true.
    self % isIsotropic  = .false.

  end subroutine kill

end module bbSurfaceSource_class
