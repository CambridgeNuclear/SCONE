module particle_class

  use numPrecision
  use universalVariables
  use genericProcedures
  use coord_class,       only : coordList
  use RNG_class,         only : RNG
  use nuclearData_inter, only : nuclearData

  implicit none
  private

  !!
  !! Particle compressed for storage
  !!
  type, public :: phaseCoord
    real(defReal)              :: wgt  = 0.0     ! Particle weight
    real(defReal),dimension(3) :: r    = 0.0     ! Global position
    real(defReal),dimension(3) :: dir  = 0.0     ! Global direction
    real(defReal)              :: E    = 0.0     ! Energy
    integer(shortInt)          :: G    = 0       ! Energy group
    logical(defBool)           :: isMG = .false. ! Is neutron multi-group
    real(defReal)              :: time = 0.0     ! Particle time position
  contains
    generic    :: assignment(=) => phaseCoord_fromParticle
    procedure  :: display       => display_phaseCoord

    procedure, private :: phaseCoord_fromParticle

  end type


  type, public :: particle
    type(coordList)            :: coords
    real(defReal)              :: E         ! Particle Energy
    integer(shortInt)          :: G         ! Particle Energy Group
    real(defReal)              :: w         ! Particle Weight
    real(defReal)              :: w0        ! Particle initial weight (for implicit, variance reduction...)

    class(RNG), pointer        :: pRNG      ! Pointer to RNG associated with the particle
    class(nuclearData),pointer :: xsData => null() ! Pointer to nuclear data

    logical(defBool)           :: isDead
    logical(defBool)           :: isMG
    real(defReal)              :: time      ! Particle time point
    real(defReal)              :: timeMax   ! Maximum neutron time before cut-off
    integer(shortInt)          :: fate = 0  ! Neutron's fate after being subjected to an operator     
    
 
  contains
     ! Build procedures
    generic              :: build => buildCE, buildMG
    generic              :: assignment(=) => particle_fromPhaseCoord

    ! Inquiry about coordinates
    procedure            :: rLocal
    procedure            :: rGlobal
    procedure            :: dirLocal
    procedure            :: dirGlobal
    procedure            :: nesting
    procedure            :: getCellIdx
    !procedure            :: getLatIdx
    procedure            :: getUniIdx
    !procedure            :: regionID
    procedure            :: matIdx

    ! Operations on coordinates
    procedure            :: moveGlobal
    procedure            :: moveLocal
    procedure            :: rotate
    procedure            :: teleport
    procedure            :: point
    procedure            :: takeAboveGeom

    ! Debug procedures
    procedure            :: display => display_particle

    !! Private - Implementation specific procedures
    procedure,private    :: buildCE
    procedure,private    :: buildMG
    procedure, private   :: particle_fromPhaseCoord

  end type particle

contains

  !!
  !! Return the position either at the deepest nested level or a specified level
  !!
  function rLocal(self,n)result(r)
    class(particle), intent(in)             :: self
    integer(shortInt), intent(in), optional :: n
    real(defReal), dimension(3)             :: r
    integer(shortInt)                       :: nMax

    if (present(n)) then
      r = self % coords % lvl(n) % r

    else
      nMax = self % coords % nesting
      r    = self % coords % lvl(nMax) % r

    end if

  end function rLocal

  !!
  !! Return the position at the highest level
  !!
  function rGlobal(self)result(r)
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
    integer(shortInt)           :: nMax

    if (present(n)) then
      dir = self % coords % lvl(n) % dir

    else
      nMax = self % coords % nesting
      dir  = self % coords % lvl(nMax) % dir

    end if

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
  !!
  function getCellIdx(self,n) result(idx)
    class(particle), intent(in)    :: self
    integer(shortInt), intent(in)  :: n
    integer(shortInt)              :: idx

    idx = self % coords % lvl(n) % cellIdx

  end function getCellIdx

!  !!
!  !! Return lattice index at a given nesting level n
!  !!
!  function getLatIdx(self,n) result(idx)
!    class(particle), intent(in)    :: self
!    integer(shortInt), intent(in)  :: n
!    integer(shortInt)              :: idx
!
!    idx = self % coords % lvl(n) % latIdx
!
!  end function getLatIdx

  !!
  !! Return universe index at a given nesting level n
  !!
  function getUniIdx(self,n) result(idx)
    class(particle), intent(in)    :: self
    integer(shortInt), intent(in)  :: n
    integer(shortInt)              :: idx

    idx = self % coords % lvl(n) % uniIdx

  end function getUniIdx

!  !!
!  !! Return unique cell id (regionID) at the lowest coordinate level
!  !!
!  function regionID(self) result(ID)
!    class(particle), intent(in) :: self
!    integer(shortInt)           :: ID
!
!    ID = self % coords % regionID
!
!  end function regionID

  !!
  !! Return current material index
  !!
  function matIdx(self) result(Idx)
    class(particle), intent(in) :: self
    integer(shortInt)           :: Idx

    Idx = self % coords % matIdx

  end function matIdx

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
  subroutine takeAboveGeom(self)
    class(particle), intent(inout) :: self

    call self % coords % takeAboveGeom()

  end subroutine takeAboveGeom

  !!
  !! Display state of a particle
  !!
  subroutine display_particle(self)
    class(particle), intent(in) :: self
    type(phaseCoord)            :: state

    state = self
    call state % display()
    print *, self % coords % matIdx

  end subroutine display_particle


  !!
  !! Copy particle into phase coordinates
  !!
  subroutine phaseCoord_fromParticle(LHS,RHS)
    class(phaseCoord), intent(out)  :: LHS
    type(particle), intent(in)      :: RHS

    LHS % wgt  = RHS % w
    LHS % r    = RHS % rGlobal()
    LHS % dir  = RHS % dirGlobal()
    LHS % E    = RHS % E
    LHS % G    = RHS % G
    LHS % isMG = RHS % isMG
    LHS % time = RHS % time

  end subroutine phaseCoord_fromParticle

  !!
  !! Copy phase coordinates into particle
  !!
  subroutine particle_fromPhaseCoord(LHS,RHS)
    class(particle), intent(inout)  :: LHS
    type(phaseCoord), intent(in)  :: RHS

    LHS % w                     = RHS % wgt
    LHS % w0                    = RHS % wgt
    call LHS % takeAboveGeom()
    LHS % coords % lvl(1) % r   = RHS % r
    LHS % coords % lvl(1) % dir = RHS % dir
    LHS % E                     = RHS % E
    LHS % G                     = RHS % G
    LHS % isMG                  = RHS % isMG
    LHS % time                  = RHS % time

  end subroutine particle_fromPhaseCoord

  !!
  !! Prints state of the phaseCoord
  !!
  subroutine display_phaseCoord(self)
    class(phaseCoord), intent(in) :: self

    print *, self % r, self % dir, self % E, self % G, self % isMG, self % wgt, self % time

  end subroutine display_phaseCoord

  !!
  !! Initialise CE particle
  !!
  subroutine buildCE(self,r,dir,E,w)
    class(particle), intent(inout)          :: self
    real(defReal),dimension(3),intent(in)   :: r, dir
    real(defReal),intent(in)                :: E, w

    call self % coords % init(r, dir)
    self % E  = E
    self % w  = w
    self % w0 = w

    self % isDead = .false.
    self % isMG   = .false.

  end subroutine buildCE

  !!
  !! Initialise MG particle
  !!
  subroutine buildMG(self,r,dir,G,w)
    class(particle), intent(inout)          :: self
    real(defReal),dimension(3),intent(in)   :: r, dir
    real(defReal),intent(in)                :: w
    integer(shortInt),intent(in)            :: G

    call self % coords % init(r, dir)
    self % G  = G
    self % w  = w
    self % w0 = w

    self % isDead = .false.
    self % isMG   = .true.

  end subroutine buildMG


end module particle_class
