module particle_class

  use numPrecision
 ! use universalVariables
  use genericProcedures
  use coord_class,    only : coordList
  use RNG_class,      only : RNG

  implicit none
  private

  type, public :: particle
    type(coordList)            :: coords
    real(defReal)              :: E         ! Particle Energy
    integer(shortInt)          :: G         ! Particle Energy Group
    real(defReal)              :: w         ! Particle Weight
    class(RNG), pointer        :: pRNG      ! Pointer to RNG associated with the particle

                               !*** Changed from currentMaterialIndex
    integer(shortInt)          :: matIdx      ! The index of the current material which the particle is traversing

    logical(defBool)           :: isDead
    logical(defBool)           :: isMG
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
    procedure            :: localDir
    procedure            :: globalDir
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
    class(particle), intent(out)            :: self
    real(defReal),dimension(3),intent(in)   :: r, dir
    real(defReal),intent(in)                :: E, w

    call self % coords % init(r, dir)
    self % E = E
    self % w = w

    self % isDead = .false.
    self % isMG = .false.
  end subroutine

  subroutine buildMG(self,r,dir,G,w)
    class(particle), intent(out)            :: self
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
    class(particle), intent(in)             :: self
    integer(shortInt), intent(in), optional :: n
    real(defReal), dimension(3)             :: r
    integer(shortInt)                       :: nMax

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
  function localDir(self,n)result(dir)
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
  end function localDir

  !!
  !! Return the direction at the global level
  !!
  function globalDir(self)result(dir)
    class(particle), intent(in) :: self
    real(defReal), dimension(3) :: dir
    dir = self % coords % lvl(1) % dir
  end function globalDir

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
  !! Rotate particle by an angle
  !!
  subroutine rotate(self,mu,phi)
    class(particle), intent(inout) :: self
    real(defReal), intent(in)      :: mu
    real(defReal), intent(in)      :: phi

    call self % coords % rotate(mu,phi)

  end subroutine rotate


  !!
  !! Place particle at an arbitrary point in the geometry in global co-ordinates
  !!
  subroutine teleport(self, r)
    class(particle), intent(inout)         :: self
    real(defReal),dimension(3), intent(in) :: r
    call self % coords % assignPosition(r)
  end subroutine teleport

  !!
  !! Point particle in an arbitrary direction in global co-ordinates
  !!
  subroutine point(self, dir)
    class(particle), intent(inout)         :: self
    real(defReal), dimension(3),intent(in) :: dir

    call self % coords % assignDirection(dir)

  end subroutine point

end module particle_class
