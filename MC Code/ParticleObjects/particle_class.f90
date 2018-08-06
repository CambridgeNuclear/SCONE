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

    class(RNG), pointer        :: pRNG      ! Pointer to RNG associated with the particle
    class(nuclearData),pointer :: xsData => null() ! Pointer to nuclear data

    logical(defBool)           :: isDead
    logical(defBool)           :: isMG
 real(defReal)              :: time      ! Particle time point
    real(defReal)              :: timeMax   ! Maximum neutron time before cut-off
    integer(shortInt)          :: fate = 0  ! Neutron's fate after being subjected to an operator     
    
 
  contains
    !! Public Interface
    generic              :: build => buildCE, buildMG
    procedure            :: makeMG
    procedure            :: moveGlobal
    procedure            :: moveLocal
    procedure            :: rotate
    procedure            :: teleport
    procedure            :: point
    procedure            :: rLocal
    procedure            :: rGlobal
    procedure            :: dirLocal
    procedure            :: dirGlobal
    procedure            :: nesting
    procedure            :: resetNesting
    procedure            :: getCellIdx
    procedure            :: getLatIdx
    procedure            :: getUniIdx
    procedure            :: updateLocation
    ! Debug procedures
    procedure            :: display => display_particle

    !! Private - Implementation specific procedures
    procedure,private    :: buildCE
    procedure,private    :: buildMG
  end type particle

contains

  subroutine makeMG(self,G)
    class(particle), intent(inout)      :: self
    integer(shortInt), intent(in)  :: G
    self % G = G
    self % isMG = .true.

  end subroutine makeMG

  subroutine buildCE(self,r,dir,E,w)
    class(particle), intent(inout)          :: self
    real(defReal),dimension(3),intent(in)   :: r, dir
    real(defReal),intent(in)                :: E, w

    call self % coords % init(r, dir)
    self % E = E
    self % w = w

    self % isDead = .false.
    self % isMG = .false.
  end subroutine

  subroutine buildMG(self,r,dir,G,w)
    class(particle), intent(inout)          :: self
    real(defReal),dimension(3),intent(in)   :: r, dir
    real(defReal),intent(in)                :: w
    integer(shortInt),intent(in)            :: G

    call self % coords % init(r, dir)
    self % G = G
    self % w = w

    self % isDead = .false.
    self % isMG = .true.
  end subroutine

  !!
  !! Return the position either at the deepest nested level or a specified level
  !!
  function rLocal(self,n)result(r)
    class(particle), intent(in) :: self
    integer(shortInt), intent(in), optional :: n
    real(defReal), dimension(3) :: r
    integer(shortInt) :: nMax
    if (present(n)) then
      r = self % coords % lvl(n) % r
    else
      nMax = self % coords % nesting
      r = self % coords % lvl(nMax) % r
    end if
  end function rLocal

  !!
  !! Return the position at the global level
  !!
  function rGlobal(self)result(r)
    class(particle), intent(in) :: self
    real(defReal), dimension(3) :: r
    r = self % coords % lvl(1) % r
  end function rGlobal

  !!
  !! Return the direction either at the deepest nested level or a specified level
  !!
  function dirLocal(self,n)result(dir)
    class(particle), intent(in) :: self
    integer(shortInt), optional :: n
    real(defReal), dimension(3) :: dir
    integer(shortInt) :: nMax
    if (present(n)) then
      dir = self % coords % lvl(n) % dir
    else
      nMax = self % coords % nesting
      dir = self % coords % lvl(nMax) % dir
    end if
  end function dirLocal

  !!
  !! Return the direction at the global level
  !!
  function dirGlobal(self)result(dir)
    class(particle), intent(in) :: self
    real(defReal), dimension(3) :: dir
    dir = self % coords % lvl(1) % dir
  end function dirGlobal

  !!
  !! Return the nestedness of the particle
  !!
  function nesting(self) result(n)
    class(particle), intent(in) :: self
    integer(shortInt)           :: n
    n = self % coords % nesting
  end function nesting

  !!
  !! Resets the particle's nesting level
  !!
  subroutine resetNesting(self)
    class(particle), intent(inout) :: self
    call self % coords % resetNesting()
  end subroutine resetNesting

  !!
  !! Get cell index for a given nested level
  !!
  function getCellIdx(self,i) result(idx)
    class(particle), intent(in)    :: self
    integer(shortInt), intent(in)  :: i
    integer(shortInt)              :: idx
    idx = self % coords % lvl(i) % cellIdx
  end function getCellIdx

  !!
  !! Get lattice index for a given nested level
  !!
  function getLatIdx(self,i) result(idx)
    class(particle), intent(in)    :: self
    integer(shortInt), intent(in)  :: i
    integer(shortInt)              :: idx
    idx = self % coords % lvl(i) % latIdx
  end function getLatIdx

  !!
  !! Get universe index for a given nested level
  !!
  function getUniIdx(self,i) result(idx)
    class(particle), intent(in)    :: self
    integer(shortInt), intent(in)  :: i
    integer(shortInt)              :: idx
    idx = self % coords % lvl(i) % uniIdx
  end function getUniIdx

  !!
  !! Set the particle's present unique region ID and material index
  !! Passes info from coordList
  !!
  subroutine updateLocation(self)
    class(particle) :: self
    self % regionID = self % coords % regionID
    self % matIdx = self % coords % matIdx
  end subroutine

  !!
  !! Move the particle in global co-ordinates only, resetting its nesting
  !!
  subroutine moveGlobal(self,distance)
    class(particle), intent(inout) :: self
    real(defReal), intent(in) :: distance
    call self % coords % moveGlobal(distance)
  end subroutine moveGlobal

  !!
  !! Move particle in local co-ordinates down to nesting level n
  !!
  subroutine moveLocal(self,distance,n)
    class(particle), intent(inout) :: self
    real(defReal), intent(in) :: distance
    integer(shortInt), intent(in) :: n
    call self % coords % moveLocal(distance,n)
  end subroutine moveLocal

  !!
  !! Place particle at an arbitrary point in the geometry in global co-ordinates
  !!
  subroutine teleport(self, r)
    class(particle), intent(inout) :: self
    real(defReal), dimension(3), intent(in) :: r
    call self % coords % assignPosition(r)
  end subroutine teleport

  !!
  !! Point particle in an arbitrary direction in global co-ordinates
  !!
  subroutine point(self, dir)
    class(particle), intent(inout) :: self
    real(defReal), dimension(3), intent(in) :: dir
    call self % coords % assignDirection(dir)
  end subroutine point

  !!
  !! Rotate particle by an angle
  !!
  subroutine rotate(self,mu,phi)
    class(particle), intent(inout) :: self
    real(defReal), intent(in)      :: mu
    real(defReal), intent(in)      :: phi

    call self % coords % rotate(mu,phi)

  end subroutine rotate


  !!
  !! Display state of a particle
  !!
  subroutine display_particle(self)
    class(particle), intent(in) :: self
    type(phaseCoord)            :: state

    state = self
    call state % display()

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
  !! Prints state of the phaseCoord
  !!
  subroutine display_phaseCoord(self)
    class(phaseCoord), intent(in) :: self

    print *, self % r, self % dir, self % E, self % G, self % isMG, self % wgt

  end subroutine display_phaseCoord

end module particle_class
