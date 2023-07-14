module pointSource_class

  use numPrecision
  use universalVariables,      only : OUTSIDE_MAT
  use genericProcedures,       only : fatalError, numToChar
  use particle_class,          only : particleState, P_NEUTRON, P_PHOTON
  use dictionary_class,        only : dictionary
  use configSource_inter,      only : configSource, kill_super => kill
  use geometry_inter,          only : geometry
  use RNG_class,               only : RNG
  use nuclearDataReg_mod,      only : ndReg_getNeutronMG => getNeutronMG
  use mgNeutronDatabase_inter, only : mgNeutronDatabase

  implicit none
  private

  !!
  !! Class describing point-like particle sources
  !!
  !! Generates a mono-directional or isotropic particle
  !! source from a single point in space, with particles of a single type
  !!
  !! MG particles may be either monoenergetic or from a distribution
  !!
  !! TODO: allow for a continuous energy distribution
  !!
  !! Private members:
  !!   r            -> source position
  !!   dir          -> optional source direction
  !!   E            -> source energy
  !!   G            -> source group
  !!   probG        -> vector of probabilities of producing particle in group G
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
  !!       type pointSource;
  !!       r (0.0 0.0 0.0);
  !!       particle neutron;
  !!       #E 14.1;              #
  !!       #G 3;                 #
  !!       #probG (0.1 0.2 0.7); #
  !!       #dir (2.0 1.0 0.0);   #
  !!      }
  !!
  type, public,extends(configSource) :: pointSource
    private
    real(defReal),dimension(3)                :: r   = ZERO
    real(defReal),dimension(3)                :: dir = ZERO
    real(defReal)                             :: E   = ZERO
    integer(shortInt)                         :: G   = 0
    real(defReal), dimension(:), allocatable  :: probG
    integer(shortInt)                         :: particleType = P_NEUTRON
    logical(defBool)                          :: isMG         = .false.
    logical(defBool)                          :: isIsotropic  = .false.
  contains
    procedure :: init
    procedure :: sampleType
    procedure :: samplePosition
    procedure :: sampleEnergy
    procedure :: sampleEnergyAngle
    procedure :: kill
  end type pointSource

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
  !!   - error if both monoenergetic and a distribution
  !!   - error if neither energy type is specified
  !!
  subroutine init(self, dict, geom)
    class(pointSource), intent(inout)      :: self
    class(dictionary), intent(in)          :: dict
    class(geometry), pointer, intent(in)   :: geom
    character(30)                          :: type
    integer(shortInt)                      :: matIdx, uniqueID
    logical(defBool)                       :: isCE, isMG, isMono, isDist
    real(defReal),dimension(:),allocatable :: temp
    class(mgNeutronDatabase), pointer      :: nucData
    character(100), parameter :: Here = 'init (pointSource_class.f90)'

    ! Provide geometry info to source
    self % geom => geom

    ! Identify which particle is used in the source
    ! Presently limited to neutron and photon
    call dict % getOrDefault(type, 'particle' ,'neutron')
    select case(type)
      case('neutron')
        self % particleType = P_NEUTRON

      case('photon')
        self % particleType = P_PHOTON

      case default
        call fatalError(Here, 'Unrecognised particle type')

    end select

    ! Get position and check it's inside geometry
    call dict % get(temp, 'r')
    if (size(self % r) /= 3) then
      call fatalError(Here, 'Source position must have three components')
    end if
    self % r = temp

    call self % geom % whatIsAt(matIdx, uniqueID, self % r)
    if (matIdx == OUTSIDE_MAT) then
      call fatalError(Here, 'Source has been placed outside geometry')
    end if

    ! Get beam direction and normalise - otherwise, assume isotropic
    self % isIsotropic = .not. dict % isPresent('dir')
    if (.not. self % isIsotropic) then

      call dict % get(temp, 'dir')
      if (size(temp) /= 3) then
        call fatalError(Here, 'Source direction must have three components')
      end if
      self % dir = temp
      self % dir = self % dir / norm2(self % dir)
    end if

    ! Get particle energy/group
    isCE = dict % isPresent('E')
    isMG = (dict % isPresent('G') .or. dict % isPresent('probG'))
    if (isCE .and. isMG) then
      call fatalError(Here, 'Source may be either continuous energy or MG, not both')

    elseif (isCE) then
      call dict % get(self % E, 'E')

    elseif (isMG) then
      isDist = dict % isPresent('probG')
      isMono = dict % isPresent('G')
      if (isDist .and. isMono) then
        call fatalError(Here, 'Source may be either monoenergetic or a distribution, not both')
      end if

      if (isDist) then
        call dict % get(temp, 'probG')
        nucData => ndReg_getNeutronMG()
        if(.not.associated(nucData)) call fatalError(Here, 'Failed to retrieve Nuclear Database')

        if (size(temp) /= nucData % nG) then
          call fatalError(Here, 'Source energy group distribution must have '//numToChar(nucData % nG)//' groups')
        end if

        allocate(self % probG(nucData % nG))
        self % probG = temp
        self % probG = self % probG / sum(self % probG)
      else
        call dict % get(self % G, 'G')
      end if

      self % isMG = .true.

    else
      call fatalError(Here, 'Must specify source energy, either Energy(E) or Group distribution (probG)')
    end if

  end subroutine init

  !!
  !! Provide particle type
  !!
  !! See configSource_inter for details.
  !!
  subroutine sampleType(self, p, rand)
    class(pointSource), intent(inout)   :: self
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
    class(pointSource), intent(inout)   :: self
    class(particleState), intent(inout) :: p
    class(RNG), intent(inout)           :: rand

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
    class(pointSource), intent(inout)   :: self
    class(particleState), intent(inout) :: p
    class(RNG), intent(inout)           :: rand
    real(defReal)                       :: r, phi, theta

    if (self % isIsotropic) then
      r = rand % get()
      phi = TWO_PI * r
      r = rand % get()
      theta = acos(1 - TWO * r)
      p % dir = [cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta)]
    else
      p % dir = self % dir
    end if

  end subroutine sampleEnergyAngle

  !!
  !! Provide particle energy
  !!
  !! See configSource_inter for details.
  !!
  subroutine sampleEnergy(self, p, rand)
    class(pointSource), intent(inout)   :: self
    class(particleState), intent(inout) :: p
    class(RNG), intent(inout)           :: rand
    real(defReal)                       :: r
    integer(shortInt)                   :: g

    if (self % isMG) then
      ! Sample from distribution
      if (allocated(self % probG)) then
        !! TODO: Replace this single, common procedure to sample discrete distributions
        r = rand % get()
        do g = 1, size(self % probG)
          r = r - self % probG(g)
          if (r < ZERO) exit
        end do
        p % G = g
        p % isMG = .true.
      ! Monoenergetic
      else
        p % G = self % G
      end if
    else
      p % E = self % E
      p % isMG = .false.
    end if

  end subroutine sampleEnergy

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(pointSource), intent(inout) :: self

    ! Kill superclass
    call kill_super(self)

    ! Kill local components
    self % r   = ZERO
    self % dir = ZERO
    self % E   = ZERO
    self % G   = 0
    if (allocated(self % probG)) deallocate(self % probG)
    self % particleType = P_NEUTRON
    self % isMG         = .false.
    self % isIsotropic  = .false.

  end subroutine kill

end module pointSource_class
