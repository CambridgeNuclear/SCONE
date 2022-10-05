module surfaceSource_class

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
  !!   - Requires deltat and nParticles in input file to be the same as specified elsewhere
  !!     in file, can change to not require these inputs with some more thought
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
  !!       type surfaceSource;
  !!       shape circle    ! circle or square;
  !!       size 5;         ! radius(circle) or side length(square)
  !!       axis x;         ! axis normal to planar shape
  !!       pos 0;          ! distance along axis to place plane
  !!       T 1;            ! temperature of source boundary
  !!       nParticles 100; ! Number of particles emitted per time step, for now has to be
  !!                         the same as IMC source if used in IMC calculation
  !!       particle photon;
  !!       #dir 1; #       ! Positive or negative to indicate direction along axis
  !!                         If 0 then emit in both directions
  !!       deltat 1;       ! Currently needed to be the same as IMC time step size
  !!      }
  !!
  type, public,extends(configSource) :: surfaceSource
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
    procedure :: kill
  end type surfaceSource

contains

  !!
  !! Initialise from dictionary
  !!
  !! See source_inter for details
  !!
  !! Errors:
  !!   - error if an unrecognised particle type is provided
  !!   - error if source is not inside geometry
  !!   - error if either direction or position have more than 3 components
  !!   - error if both CE and MG is specified
  !!   - error if neither energy type is specified
  !!
  subroutine init(self, dict, geom)
    class(surfaceSource), intent(inout)    :: self
    class(dictionary), intent(in)          :: dict
    class(geometry), pointer, intent(in)   :: geom
    character(30)                          :: type, tempName
    integer(shortInt)                      :: matIdx, uniqueID
    logical(defBool)                       :: isCE, isMG
    real(defReal) :: temp !,dimension(:),allocatable :: temp
    character(100), parameter :: Here = 'init (surfaceSource_class.f90)'

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
        call fatalError(Here, 'Unrecognised axis, may onlt be x, y or z')
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
    call dict % get(self % deltat, 'deltat')

  end subroutine init

  subroutine append(self, dungeon, N, rand)
    class(surfaceSource), intent(inout)  :: self
    type(particleDungeon), intent(inout) :: dungeon
    integer(shortInt), intent(in)        :: N
    class(RNG), intent(inout)            :: rand
    integer(shortInt)                    :: i
    character(100), parameter            :: Here = 'append (surfaceSource_class.f90)'

    self % N = N

    ! Generate n particles to populate dungeon
    do i = 1, N
      call dungeon % detain(self % sampleParticle(rand))
    end do

  end subroutine append

  !!
  !! Provide particle type
  !!
  !! See configSource_inter for details.
  !!
  subroutine sampleType(self, p, rand)
    class(surfaceSource), intent(inout)   :: self
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
    class(surfaceSource), intent(inout)   :: self
    class(particleState), intent(inout) :: p
    class(RNG), intent(inout)           :: rand
    real(defReal), dimension(3)         :: prevPos
    real(defReal)                       :: r1, r2, rad, theta

    if ( self % planeShape == 0 ) then   ! Square

      prevPos = self % r

      ! Set new x, y and z coords
      self % r(1) = rand % get() * self % surfSize/2
      self % r(2) = rand % get() * self % surfSize/2
      self % r(3) = rand % get() * self % surfSize/2
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
    class(surfaceSource), intent(inout)   :: self
    class(particleState), intent(inout) :: p
    class(RNG), intent(inout)           :: rand
    real(defReal)                       :: r, phi, theta

    r = rand % get()
    phi = TWO_PI * r
    r = rand % get()
    theta = acos(1 - TWO * r)
    p % dir = [cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta)]

    ! If dir not equal to zero, adjust so that particles are travelling in correct direction
    if (self % dir /= 0) then
      p % dir(self % axis) = abs(p % dir(self % axis)) * self % dir 
    end if

  end subroutine sampleEnergyAngle

  !!
  !! Provide particle energy
  !!
  !! Sampled as a black body surface, see "Four Decades of Implicit Monte Carlo",
  !!  Allan B Wollaber, p.24-25
  !!
  !! See configSource_inter for details.
  !!
  subroutine sampleEnergy(self, p, rand)
    class(surfaceSource), intent(inout) :: self
    class(particleState), intent(inout) :: p
    class(RNG), intent(inout)           :: rand
    real(defReal)                       :: num

    num = radiationConstant * lightSpeed * self % deltat * self % T**4 * self % area
    p % wgt = num / (4 * self % N)

    ! If dir = 0 then emit in both directions => double total energy
    if (self % dir == 0) p % wgt = 2*p % wgt

    p % isMG = .true.
    p % G = 1

  end subroutine sampleEnergy

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(surfaceSource), intent(inout) :: self

    ! Kill superclass
    call kill_super(self)

    ! Kill local components
    self % r   = ZERO
    self % dir = ZERO
    self % particleType = P_PHOTON
    self % isMG         = .true.
    self % isIsotropic  = .false.

  end subroutine kill

end module surfaceSource_class
