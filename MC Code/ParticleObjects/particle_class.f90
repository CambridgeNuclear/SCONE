module particle_class

  use numPrecision
 ! use universalVariables
  use genericProcedures, only : rotateVector
  !use cell_class
  !use surface_class
  use coord_class,    only : coord, coordList
  use RNG_class,      only : RNG


  implicit none
  private

  type, public :: particle
    type(coordList)            :: coords
    real(defReal)              :: E         ! Particle Energy
    integer(shortInt)          :: G         ! Particle Energy Group
    real(defReal)              :: w         ! Particle Weight
    real(defReal),dimension(3) :: dir       ! *** Debug dummy direction

   ! type(cell_ptr)             :: currentCell               ! The current cell which the particle occupies
   ! class(surface), pointer    :: currentSurface => null()  ! The current surface on which the particle is sat
    integer(shortInt)          :: matIdx      ! The index of the current material which the particle is traversing

    class(RNG), pointer        :: pRNG     ! Pointer to RNG associated with the particle

    logical(defBool)           :: isDead
    logical(defBool)           :: isMG
  contains
    !! Public Interface
   ! generic              :: build => buildCE, buildMG
    procedure             :: makeMG
    procedure             :: rotate
    procedure             :: getDirection
    procedure             :: point
    !procedure            :: turnGlobal
   ! procedure            :: moveGlobal
   ! procedure            :: moveLocal
    !! Private - Implementation specific procedures
   ! procedure,private    :: buildCE
   ! procedure,private    :: buildMG
  end type particle

contains

  function getDirection(self) result(dir)
    class(particle), intent(in)  :: self
    real(defReal), dimension(3)  :: dir

    dir = self % dir

  end function getDirection

  subroutine rotate(self,mu,phi)
    class(particle), intent(inout) :: self
    real(defReal), intent(in)      :: mu
    real(defReal), intent(in)      :: phi

    self % dir = rotateVector(self % dir,mu,phi)

  end subroutine rotate

  subroutine point(self,dir)
    class(particle), intent(inout)        :: self
    real(defReal),dimension(3),intent(in) :: dir
    self % dir = dir

  end subroutine point



  subroutine makeMG(self,G)
    class(particle), intent(inout)      :: self
    integer(shortInt), intent(in)  :: G
    self % G = G
    self % isMG = .true.

  end subroutine makeMG
!
!  subroutine buildCE(self,r,dir,E,w)
!    class(particle), intent(out)            :: self
!    real(defReal),dimension(3),intent(in)   :: r, dir
!    real(defReal),intent(in)                :: E, w
!
!    call self % coords % init(r, dir)
!    self % E = E
!    self % w = w
!
!    self % isDead = .false.
!    self % isMG = .false.
!  end subroutine
!
!  subroutine buildMG(self,r,dir,G,w)
!    class(particle), intent(out)            :: self
!    real(defReal),dimension(3),intent(in)   :: r, dir
!    real(defReal),intent(in)                :: w
!    integer(shortInt),intent(in)            :: G
!
!    call self % coords % init(r, dir)
!    self % G = G
!    self % w = w
!
!    self % isDead = .false.
!    self % isMG = .true.
!  end subroutine

  !!
  !! Change the direction of a particle by provided cosines of polar and azimuthal deflection angle
  !!
  !subroutine turnGlobal(self


  !!
  !! Move the particle in global co-ordinates only, resetting its nesting
  !!
!  subroutine moveGlobal(self,distance)
!    implicit none
!    class(particle), intent(inout) :: self
!    real(defReal), intent(in) :: distance
!
!    call self % coords % resetNesting()
!    self % coords % lvl(1) % r = self % coords % lvl(1) % r + &
!                                 self % coords % lvl(1) % dir * distance
!  end subroutine moveGlobal

  !!
  !! Move particle in local co-ordinates down to nesting level n
  !!
!  subroutine moveLocal(self,distance,n)
!    implicit none
!    class(particle), intent(inout) :: self
!    real(defReal), intent(in) :: distance
!    integer(shortInt), intent(in) :: n
!    integer(shortInt) :: i
!
!    call self % coords % resetNesting(n)
!    do i = 1,n
!      self % coords % lvl(i) % r = self % coords % lvl(i) % r + &
!                                 self % coords % lvl(i) % dir * distance
!    end do
!  end subroutine moveLocal

end module particle_class
