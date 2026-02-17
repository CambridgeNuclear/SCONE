module particle_class

  use numPrecision
  use universalVariables
  use genericProcedures
  use coord_class,       only : coordList
  use RNG_class,         only : RNG
  use errors_mod,        only : fatalError
  use tallyCodes

  implicit none
  private

  !!
  !! Particle types parameters
  !!
  integer(shortInt), parameter,public :: P_NEUTRON   = 1, &
                                         P_PHOTON    = 2, &
                                         P_PRECURSOR = 3

  !!
  !! Public particle type procedures
  !!
  public :: verifyType
  public :: printType

  !!
  !! Particle compressed for storage
  !!
  !! particleStateData contains only the properties (useful for MPI); it is extended into
  !! particleState which includes the procedures too and is used for tallying.
  !!
  !! Public Members:
  !!   wgt      -> Weight of the particle
  !!   r        -> Global Position of the particle [cm]
  !!   dir      -> Direction vector of the particle (normalised to 1.0)
  !!   E        -> Energy of the particle [MeV]
  !!   G        -> Energy Group of the particle
  !!   isMG     -> True if particle uses MG data
  !!   type     -> Physical Type of the particle (NEUTRON, PHOTON etc.)
  !!   time     -> Position in time of the particle [s]
  !!   matIdx   -> material Index in which particle is present
  !!   cellIdx  -> Cell Index at the lowest level in which particle is present
  !!   uniqueID -> Unique ID of the cell at the lowest level in which particle is present
  !!   collisionN -> Number of collisions the particle went through
  !!   broodID  -> ID of the source particle. It is used to indicate the primogenitor of the particles
  !!               in the particleDungeon so they can be sorted, which is necessary for reproducibility
  !!               with OpenMP
  !!
  type, public :: particleStateData
    real(defReal)              :: wgt  = ZERO       ! Particle weight
    real(defReal),dimension(3) :: r    = ZERO       ! Global position
    real(defReal),dimension(3) :: dir  = ZERO       ! Global direction
    real(defReal)              :: E    = ZERO       ! Energy
    integer(shortInt)          :: G    = 0          ! Energy group
    logical(defBool)           :: isMG = .false.    ! Is neutron multi-group
    integer(shortInt)          :: type = P_NEUTRON  ! Particle physical type
    real(defReal)              :: time = ZERO       ! Particle time position
    real(defReal)              :: lambda = INF      ! Precursor decay constant
    integer(shortInt)          :: matIdx   = -1     ! Material index where particle is
    integer(shortInt)          :: cellIdx  = -1     ! Cell idx at the lowest coord level
    integer(shortInt)          :: uniqueID = -1     ! Unique id at the lowest coord level
    integer(shortInt)          :: collisionN = 0    ! Number of collisions
    integer(shortInt)          :: broodID = 0       ! ID of the source particle
  end type particleStateData

  !!
  !! Extension of particleStateData, which includes procedure
  !!
  !! Interface:
  !!   assignemnt(=)  -> Build particleState from particle
  !!   operator(.eq.) -> Return True if particle are exactly the same
  !!   display        -> Print debug information about the state to the console
  !!
  type, public, extends(particleStateData) :: particleState
  contains
    generic    :: assignment(=)  => fromParticle
    generic    :: assignment(=)  => fromParticleStateData
    generic    :: operator(.eq.) => equal_particleState
    procedure  :: display        => display_particleState
    procedure  :: kill           => kill_particleState
    procedure  :: fromParticle   => particleState_fromParticle
    procedure  :: fromParticleStateData => particleState_fromParticleStateData

    ! Private procedures
    procedure, private :: equal_particleState

  end type particleState

  !!
  !! This type represents particle
  !!
  !! In current form it was designed to support neutron and other neutral particles
  !! By extension(inheritance) support for photons or charged particles could be introduced
  !!
  !! In addition to representing a particle by physical parameters, it containes additional
  !! data like RNG pointer. This is to enable access to objects associated with a particlular
  !! particle history from anywhere where a particle was passed. This it works kind of like a
  !! global scope for a specific particle. Such setup should enable particle to control how it
  !! is beeing processed during Monte Carlo transport.
  !!
  !! Particle also contains snapshots (states) at privious points in its history. This is used
  !! for tallying. It is necessary to compare current state with the earlier state to score e.g.
  !! collision probabilities.
  !!
  type, public :: particle
    ! Particle phase space data
    type(coordList)            :: coords
    real(defReal)              :: E         ! Particle Energy
    integer(shortInt)          :: G         ! Particle Energy Group
    real(defReal)              :: w         ! Particle Weight
    real(defReal)              :: time      ! Particle time point

    ! Precursor particle data
    real(defReal)              :: lambda = INF     ! Precursor decay constant
    
    ! Particle flags
    real(defReal)              :: w0             ! Particle initial weight (for implicit, variance reduction...)
    logical(defBool)           :: isDead
    logical(defBool)           :: isMG
    real(defReal)              :: timeMax = -INF ! Maximum neutron time before cut-off
    integer(shortInt)          :: fate = NO_FATE ! Neutron's fate after being subjected to an operator
    integer(shortInt)          :: type           ! Particle type
    integer(shortInt)          :: collisionN = 0 ! Index of the number of collisions the particle went through
    integer(shortInt)          :: broodID = 0    ! ID of the brood (source particle number)

    ! Particle processing information
    class(RNG), pointer        :: pRNG  => null()   ! Pointer to RNG associated with the particle
    real(defReal)              :: k_eff             ! Value of default keff for implicit source generation
    real(defReal)              :: alpha = ZERO      ! Value of alpha for time absorption/production
    real(defReal)              :: eta = ONE         ! Value for stabilising alpha iteration
    integer(shortInt)          :: geomIdx           ! Index of the geometry used by the particle
    integer(shortInt)          :: splitCount = 0    ! Counter of number of splits

    ! Archived snapshots of previous states
    type(particleState)        :: preHistory
    type(particleState)        :: preTransition
    type(particleState)        :: prePath
    type(particleState)        :: preCollision

  contains
     ! Build procedures
    generic              :: build => buildCE, buildMG
    generic              :: assignment(=) => particle_fromParticleState

    ! Enquiry about coordinates
    procedure                  :: rLocal
    procedure                  :: rGlobal
    procedure                  :: dirLocal
    procedure                  :: dirGlobal
    procedure                  :: nesting
    procedure                  :: getCellIdx
    procedure                  :: getUniIdx
    procedure                  :: matIdx
    procedure, non_overridable :: getType

    ! Enquiry about physical state
    procedure :: getSpeed
    procedure :: getAlphaAbsorption

    ! Precursor procedures
    procedure :: isPrecursor
    procedure :: emitDelayedNeutron
    procedure :: forcedPrecursorDecay
    procedure :: expectedDelayedWgt

    ! Operations on coordinates
    procedure :: moveGlobal
    procedure :: moveLocal
    procedure :: rotate
    procedure :: teleport
    procedure :: point
    procedure :: takeAboveGeom
    procedure :: setMatIdx

    ! Save particle state information
    procedure, non_overridable  :: savePreHistory
    procedure, non_overridable  :: savePreTransition
    procedure, non_overridable  :: savePrePath
    procedure, non_overridable  :: savePreCollision

    ! Debug procedures
    procedure            :: display => display_particle
    procedure            :: typeToChar

    !! Private - Implementation specific procedures
    procedure,private                   :: buildCE
    procedure,private                   :: buildMG
    procedure,non_overridable,private   :: particle_fromParticleState

  end type particle

contains

!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Particle build and assignment procedures
!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Initialise CE particle
  !! Necessary arguments:
  !!   r   -> Global Position
  !!   dir -> Global direction
  !!   E   -> Energy [MeV]
  !!   w   -> particle weight
  !! Optional arguments:
  !!   t   -> particle time (default = 0.0)
  !!   type-> particle type (default = P_NEUTRON)
  !!
  pure subroutine buildCE(self, r, dir, E, w, t, type)
    class(particle), intent(inout)          :: self
    real(defReal),dimension(3),intent(in)   :: r
    real(defReal),dimension(3),intent(in)   :: dir
    real(defReal),intent(in)                :: E
    real(defReal),intent(in)                :: w
    real(defReal),optional,intent(in)       :: t
    integer(shortInt),intent(in),optional   :: type

    call self % coords % init(r, dir)
    self % E  = E
    self % w  = w
    self % w0 = w

    self % isDead = .false.
    self % isMG   = .false.

    if(present(t)) then
      self % time = t
    else
      self % time = ZERO
    end if

    if(present(type)) then
      self % type = type
    else
      self % type = P_NEUTRON
    end if

  end subroutine buildCE

  !!
  !! Initialise MG particle
  !! Necessary arguments:
  !!   r   -> Global Position
  !!   dir -> Global direction
  !!   G   -> Energy Group
  !!   w   -> particle weight
  !! Optional arguments:
  !!   t   -> particle time (default = 0.0)
  !!   type-> particle type (default = P_NEUTRON)
  !!
  subroutine buildMG(self, r, dir, G, w, t, type)
    class(particle), intent(inout)          :: self
    real(defReal),dimension(3),intent(in)   :: r
    real(defReal),dimension(3),intent(in)   :: dir
    real(defReal),intent(in)                :: w
    integer(shortInt),intent(in)            :: G
    real(defReal),intent(in),optional       :: t
    integer(shortInt),intent(in),optional   :: type

    call self % coords % init(r, dir)
    self % G  = G
    self % w  = w
    self % w0 = w

    self % isDead = .false.
    self % isMG   = .true.

    if(present(t)) then
      self % time = t
    else
      self % time = ZERO
    end if

    if(present(type)) then
      self % type = type
    else
      self % type = P_NEUTRON
    end if

  end subroutine buildMG

  !!
  !! Copy phase coordinates into particle
  !!
  pure subroutine particle_fromParticleState(LHS,RHS)
    class(particle), intent(inout)   :: LHS
    type(particleState), intent(in)  :: RHS

    LHS % w                     = RHS % wgt
    LHS % w0                    = RHS % wgt
    call LHS % takeAboveGeom()
    LHS % coords % lvl(1) % r   = RHS % r
    LHS % coords % lvl(1) % dir = RHS % dir
    LHS % E                     = RHS % E
    LHS % G                     = RHS % G
    LHS % isMG                  = RHS % isMG
    LHS % type                  = RHS % type
    LHS % time                  = RHS % time
    LHS % lambda                = RHS % lambda
    LHS % fate                  = NO_FATE
    LHS % collisionN            = RHS % collisionN
    LHS % splitCount            = 0 ! Reinitialise counter for number of splits
    LHS % broodID               = RHS % broodID

  end subroutine particle_fromParticleState

!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Particle coordinates inquiry procedures
!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Return the position either at the deepest nested level or a specified level
  !!
  function rLocal(self,n)result(r)
    class(particle), intent(in)             :: self
    integer(shortInt), intent(in), optional :: n
    real(defReal), dimension(3)             :: r
    integer(shortInt)                       :: n_loc

    if(present(n)) then
      n_loc = n
    else
      n_loc = self % coords % nesting
    end if

    r = self % coords % lvl(n_loc) % r

  end function rLocal

  !!
  !! Return the position at the highest level
  !!
  pure function rGlobal(self)result(r)
    class(particle), intent(in) :: self
    real(defReal), dimension(3) :: r

    r = self % coords % lvl(1) % r

  end function rGlobal

  !!
  !! Return the direction either at the deepest nested level or at a specified level
  !!
  function dirLocal(self,n)result(dir)
    class(particle), intent(in) :: self
    integer(shortInt), optional :: n
    real(defReal), dimension(3) :: dir
    integer(shortInt)           :: n_loc

    if(present(n)) then
      n_loc = n
    else
      n_loc = self % coords % nesting
    end if

    dir = self % coords % lvl(n_loc) % dir

  end function dirLocal

  !!
  !! Return the direction at the highest nesting level
  !!
  function dirGlobal(self)result(dir)
    class(particle), intent(in) :: self
    real(defReal), dimension(3) :: dir

    dir = self % coords % lvl(1) % dir

  end function dirGlobal

  !!
  !! Return the lowest nesting level of the particle
  !!
  function nesting(self) result(n)
    class(particle), intent(in) :: self
    integer(shortInt)           :: n

    n = self % coords % nesting

  end function nesting

  !!
  !! Return cell index at a given nesting level n
  !!  If no n is given return for lowest nesting level
  !!
  pure function getCellIdx(self,n) result(idx)
    class(particle), intent(in)             :: self
    integer(shortInt),optional, intent(in)  :: n
    integer(shortInt)                       :: idx
    integer(shortInt)                       :: n_loc

    if(present(n)) then
      n_loc = n
    else
      n_loc = self % coords % nesting
    end if

    idx = self % coords % lvl(n_loc) % cellIdx

  end function getCellIdx

  !!
  !! Return universe index at a given nesting level n
  !!
  pure function getUniIdx(self,n) result(idx)
    class(particle), intent(in)             :: self
    integer(shortInt),optional, intent(in)  :: n
    integer(shortInt)                       :: idx
    integer(shortInt)                       :: n_loc

    if(present(n)) then
      n_loc = n
    else
      n_loc = self % coords % nesting
    end if

    idx = self % coords % lvl(n_loc) % uniIdx

  end function getUniIdx

  !!
  !! Return current material index
  !!
  pure function matIdx(self) result(Idx)
    class(particle), intent(in) :: self
    integer(shortInt)           :: idx

    idx = self % coords % matIdx

  end function matIdx

  !!
  !! Return one of the particle Tpes defined in universal variables
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   P_NEUTRON_CE, P_NEUTRON_MG
  !!
  !! Errors:
  !!   None
  !!
  pure function getType(self) result(type)
    class(particle), intent(in) :: self
    integer(shortInt)           :: type

    if (self % isMG) then
      type = P_NEUTRON_MG
    else
      type = P_NEUTRON_CE
    end if

  end function getType

  !!
  !! Return the particle speed in [cm/s]
  !! neutronMass: [MeV]
  !! lightSpeed:  [cm/s]
  !!
  !! NOTE:
  !!   The speeds are computed from non-relativistic formula for massive particles.
  !!   A small error might appear in MeV range (e.g. for fusion applications).
  !!   Further there is currently no good solution for MG neutrons. Their speed
  !!   is arbitrarily set to 1.
  !!
  !! Does not provide errors if a dubious particle type is used. This is to allow
  !! the function and related functions to be pure.
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   Particle speed
  !!
  pure function getSpeed(self) result(speed)
    class(particle), intent(in) :: self
    real(defReal)               :: speed

    ! Calculates the velocity for the relevant particle [cm/s]
    if (self % type == P_PHOTON) then
      speed = lightSpeed

    elseif (self % type == P_NEUTRON) then
      ! TODO: handle MG speed
      if (self % isMG) then
        speed = ONE
      else
        speed = sqrt(TWO * self % E / neutronMass) * lightSpeed
      end if

    else
      speed = ONE

    end if

  end function getSpeed
  
  !!
  !! Return the alpha absorption cross section: XS = alpha/v
  !! Also includes Zoia's stabilisation for negative alpha values.
  !!
  pure function getAlphaAbsorption(self) result(xs)
    class(particle), intent(in) :: self
    real(defReal)               :: xs
    real(defReal)               :: alpha, eta
    
    xs = ZERO
    alpha = self % alpha
    eta = self % eta
    xs = (max(alpha, ZERO) + eta * max(-alpha, ZERO)) / self % getSpeed()

  end function getAlphaAbsorption

  !!
  !! Produce a delayed neutron from a precursor.
  !! Handles changes of type, removal of lambda.
  !!
  !! Errors:
  !!   Fatal error if producing particle is not a precursor
  !!
  subroutine emitDelayedNeutron(self)
    class(particle), intent(inout) :: self
    character(100), parameter      :: Here = 'emitDelayedNeutron (particle_class.f90)'
    
    ! Ensure decaying particle is a precursor
    if (self % type /= P_PRECURSOR) then
      call fatalError(Here, 'Only precursors (3) can decay to neutrons. Particle type is '&
              //numToChar(self % type))
    end if
    
    self % type = P_NEUTRON
    self % lambda = INF

  end subroutine emitDelayedNeutron

  !!
  !! Produces a neutron of appropriate weight given a decay time.
  !! Must be produced from a precursor.
  !!
  !! Args:
  !!   t      -> time of forced decay
  !!   deltaT -> time width across which forced decay occurs
  !!
  !! Result:
  !!   Neutron with appropriate weight at time t
  !!
  !! Errors:
  !!   Fatal error if decaying particle is not a precursor.
  !!   Fatal error if decay time is less than precursor time
  !!   Fatal error if time increment is negative
  !!
  subroutine forcedPrecursorDecay(self, t, deltaT, delayedN)
    class(particle), intent(in) :: self
    real(defReal), intent(in)   :: t
    real(defReal), intent(in)   :: deltaT
    type(particle), intent(out) :: delayedN
    real(defReal)               :: wDelay
    character(100), parameter   :: Here = 'forcedPrecursorDecay (particle_class.f90)'

    ! Ensure decaying particle is a precursor
    if (self % type /= P_PRECURSOR) then
      call fatalError(Here, 'Can only perform decay on a precursor (3). Particle type is '&
              //numToChar(self % type))
    end if
    if (t < self % time) then
      call fatalError(Here, 'Decay time must come after the precursor production time. &
             & Decay time is: '//numToChar(t)//' Production time is: '//numToChar(self % time))
    end if
    if (deltaT < ZERO) call fatalError(Here,'Time increment must be positive: '//numToChar(deltaT))

    wDelay = self % w * deltaT * self % lambda * exp(-self % lambda * (t - self % time))

    delayedN = self
    call delayedN % emitDelayedNeutron()
    delayedN % w = wDelay
    delayedN % time = t

  end subroutine forcedPrecursorDecay
  
  !!
  !! Return expected weight of delayed neutron across
  !! time interval [t1, t2]
  !!
  !! Args:
  !!   t1  -> initial time (beginning of step)
  !!   t2  -> end time (end of step)
  !!
  !! Result:
  !!   Expected delayed neutron weight
  !!
  !! Errors:
  !!   Return an error if the particle is not a precursor
  !!
  function expectedDelayedWgt(self, t1, t2) result(wgt)
    class(particle), intent(in) :: self
    real(defReal), intent(in)   :: t1
    real(defReal), intent(in)   :: t2
    real(defReal)               :: wgt
    real(defReal)               :: lam, t0
    character(100), parameter   :: Here = 'expectedDelayedWgt (particle_class.f90)'

    ! Ensure decaying particle is a precursor
    if (self % type /= P_PRECURSOR) then
      call fatalError(Here, 'Can only estimate decay weight of a precursor (3). Particle type is '&
              //numToChar(self % type))
    end if

    t0 = self % time
    
    ! Ensure sensible times used
    if (t1 >= t2) call fatalError(Here, 't1 must be less than t2. t1: '//numToChar(t1)//&
            ' t2: '//numToChar(t2))
    if (t1 < t0) call fatalError(Here, 't1 must be greater than particle time. t1: '&
            //numToChar(t1)//' Particle time: '//numToChar(t0))

    lam = self % lambda

    wgt = self % w * (exp(-lam * (t1 - t0)) - exp(-lam * (t2 - t0)))

  end function expectedDelayedWgt

  !!
  !! Returns whether the particle is a precursor
  !!
  function isPrecursor(self) result(isPrec)
    class(particle), intent(in) :: self
    logical(defBool)            :: isPrec

    isPrec = self % type == P_PRECURSOR

  end function isPrecursor

!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Particle operations on coordinates procedures
!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Move the particle above the geometry
  !! NOTE: regionID & matIdx will be reset!!!
  !!
  subroutine moveGlobal(self,distance)
    class(particle), intent(inout) :: self
    real(defReal), intent(in)      :: distance

    call self % coords % moveGlobal(distance)

  end subroutine moveGlobal

  !!
  !! Move particle in local co-ordinates down to nesting level n
  !!
  subroutine moveLocal(self,distance,n)
    class(particle), intent(inout) :: self
    real(defReal), intent(in)      :: distance
    integer(shortInt), intent(in)  :: n

    call self % coords % moveLocal(distance,n)

  end subroutine moveLocal

  !!
  !! Rotate particle
  !!  mu  -> cosine of deflection from current direction
  !!  phi -> azimuthal angle of rotation
  !!
  subroutine rotate(self,mu,phi)
    class(particle), intent(inout) :: self
    real(defReal), intent(in)      :: mu
    real(defReal), intent(in)      :: phi

    call self % coords % rotate(mu,phi)

  end subroutine rotate

  !!
  !! Place particle at an arbitrary point in above the geometry
  !!
  subroutine teleport(self, r)
    class(particle), intent(inout)          :: self
    real(defReal), dimension(3), intent(in) :: r

    call self % coords % assignPosition(r)

  end subroutine teleport

  !!
  !! Point particle in direction dir in highest nesting level
  !! Propagates new direction to lower levels
  !!
  subroutine point(self, dir)
    class(particle), intent(inout)          :: self
    real(defReal), dimension(3), intent(in) :: dir

    call self % coords % assignDirection(dir)

  end subroutine point

  !!
  !! Resets the particle's nesting level
  !!
  pure subroutine takeAboveGeom(self)
    class(particle), intent(inout) :: self

    call self % coords % takeAboveGeom()

  end subroutine takeAboveGeom

  !!
  !! Set Material index for testing purposes
  !!
  pure subroutine setMatIdx(self,matIdx)
    class(particle), intent(inout) :: self
    integer(shortInt), intent(in)  :: matIdx

    self % coords % matIdx = matIdx

  end subroutine setMatIdx

!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Particle save state procedures
!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Save state of the particle at the beginning of history
  !!
  subroutine savePreHistory(self)
    class(particle), intent(inout) :: self

    self % preHistory = self

  end subroutine savePreHistory

  !!
  !! Save state of the particle at the beginning of history
  !!
  subroutine savePreTransition(self)
    class(particle), intent(inout) :: self

    self % preTransition = self

  end subroutine savePreTransition

  !!
  !! Save state of the particle at the beginning of history
  !!
  subroutine savePrePath(self)
    class(particle), intent(inout) :: self

    self % prePath = self

  end subroutine savePrePath

  !!
  !! Save state of the particle at the beginning of history
  !!
  subroutine savePreCollision(self)
    class(particle), intent(inout) :: self

    self % preCollision = self

  end subroutine savePreCollision

!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Particle debug procedures
!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Display state of a particle
  !!
  subroutine display_particle(self)
    class(particle), intent(in) :: self
    type(particleState)         :: state

    state = self
    call state % display()
    print *, 'Material: ', self % coords % matIdx

  end subroutine display_particle

  !!
  !! Return character that describes the type of particle
  !!
  function typeToChar(self) result(c)
    class(particle), intent(in) :: self
    character(:), allocatable   :: c
    character(2)                :: eType

    if( self % isMG) then
      eType = 'MG'
    else
      eType = 'CE'
    end if

    c = eType // ' ' // trim(printType( self % type))

  end function typeToChar


!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Particle state and phaseCoord procedures
!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Copy particle into phase coordinates
  !!
  subroutine particleState_fromParticle(LHS,RHS)
    class(particleState), intent(out)  :: LHS
    class(particle), intent(in)        :: RHS

    LHS % wgt  = RHS % w
    LHS % r    = RHS % rGlobal()
    LHS % dir  = RHS % dirGlobal()
    LHS % E    = RHS % E
    LHS % G    = RHS % G
    LHS % isMG = RHS % isMG
    LHS % type = RHS % type
    LHS % time = RHS % time
    LHS % lambda = RHS % lambda

    ! Save all indexes
    LHS % matIdx   = RHS % coords % matIdx
    LHS % uniqueID = RHS % coords % uniqueId
    LHS % cellIdx  = RHS % coords % cell()
    LHS % collisionN = RHS % collisionN
    LHS % broodID    = RHS % broodID

  end subroutine particleState_fromParticle

  !!
  !! Copy particleStateData into phase coordinates
  !!
  subroutine particleState_fromParticleStateData(LHS,RHS)
    class(particleState), intent(out)    :: LHS
    class(particleStateData), intent(in) :: RHS

    LHS % wgt  = RHS % wgt
    LHS % r    = RHS % r
    LHS % dir  = RHS % dir
    LHS % E    = RHS % E
    LHS % G    = RHS % G
    LHS % isMG = RHS % isMG
    LHS % type = RHS % type
    LHS % time = RHS % time

    ! Save all indexes
    LHS % matIdx   = RHS % matIdx
    LHS % uniqueID = RHS % uniqueId
    LHS % cellIdx  = RHS % cellIdx
    LHS % collisionN = RHS % collisionN
    LHS % broodID    = RHS % broodID

  end subroutine particleState_fromParticleStateData

  !!
  !! Define equal operation on phase coordinates
  !!  Phase coords are equal if all their components are the same
  !!
  function equal_particleState(LHS,RHS) result(isEqual)
    class(particleState), intent(in) :: LHS
    class(particleState), intent(in) :: RHS
    logical(defBool)              :: isEqual

    isEqual = .true.
    isEqual = isEqual .and. LHS % wgt == RHS % wgt
    isEqual = isEqual .and. all(LHS % r   == RHS % r)
    isEqual = isEqual .and. all(LHS % dir == RHS % dir)
    isEqual = isEqual .and. LHS % time == RHS % time
    isEqual = isEqual .and. LHS % isMG .eqv. RHS % isMG
    isEqual = isEqual .and. LHS % type == RHS % type
    isEqual = isEqual .and. LHS % matIdx   == RHS % matIdx
    isEqual = isEqual .and. LHS % cellIdx  == RHS % cellIdx
    isEqual = isEqual .and. LHS % uniqueID == RHS % uniqueID
    isEqual = isEqual .and. LHS % collisionN == RHS % collisionN
    isEqual = isEqual .and. LHS % broodID    == RHS % broodID

    if( LHS % isMG ) then
      isEqual = isEqual .and. LHS % G == RHS % G
    else
      isEqual = isEqual .and. LHS % E == RHS % E
    end if
  end function equal_particleState

  !!
  !! Prints state of the phaseCoord
  !!
  subroutine display_particleState(self)
    class(particleState), intent(in) :: self

    print*, 'Position: ', self % r
    print*, 'Direction: ', self % dir
    print*, 'Energy: ', self % E
    print*, 'Group: ', self % G
    print*, 'isMG: ', self % isMG
    print*, 'Weight: ', self % wgt
    print*, 'Time: ', self % time

  end subroutine display_particleState

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill_particleState(self)
    class(particleState), intent(inout) :: self

    self % wgt  = ZERO
    self % r    = ZERO
    self % dir  = ZERO
    self % E    = ZERO
    self % G    = 0
    self % isMG = .false.
    self % type = P_NEUTRON
    self % time = ZERO
    self % matIdx   = -1
    self % cellIdx  = -1
    self % uniqueID = -1
    self % collisionN = 0
    self % broodID    = 0

  end subroutine kill_particleState


!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Misc Procedures
!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  !!
  !! Returns .true. if the integer is valid particle type code
  !! Returns .false. otherwise
  !!
  elemental function verifyType(type) result(isValid)
    integer(shortInt), intent(in) :: type
    logical(defBool)              :: isValid

    isValid = .false.

    ! Check against particles types
    isValid = isValid .or. type == P_NEUTRON
    isValid = isValid .or. type == P_PHOTON
    isValid = isValid .or. type == P_PRECURSOR

  end function verifyType

  !!
  !! Returns character with a description of particle type
  !!
  pure function printType(type) result(name)
    integer(shortInt), intent(in) :: type
    character(:),allocatable      :: name

    select case(type)
      case(P_NEUTRON)
        name = 'Neutron'

      case(P_PHOTON)
        name = 'Photon'

      case default
        name = 'INVALID PARTICLE TYPE'

    end select
  end function printType

end module particle_class
