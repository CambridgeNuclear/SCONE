module particle_class

  use numPrecision

  implicit none
  private

  type, public :: particle
    private
    real(kind=defReal),dimension(3) :: r         ! Particle Position
    real(kind=defReal),dimension(3) :: dir       ! Particle Direction
    real(kind=defReal)              :: E         ! Particle Energy
    integer(kind=shortInt)          :: G         ! Particle Energy Group
    real(kind=defReal)              :: w         ! Particle Weight

    integer(kind=shortInt)          :: CurrentCellIndex
    integer(kind=shortInt)          :: CurrentMaterialIndex

    logical(kind=defBool)           :: isDead
    logical(kind=defBool)           :: isMG
  contains
        !! Public Interface
        generic              :: build => buildCE, buildMG
        procedure            :: makeMG
        !! Private - Implementation specific procedures
        procedure,private    :: buildCE
        procedure,private    :: buildMG
  end type particle




contains

  subroutine makeMG(self,G)
    class(particle), intent(inout)      :: self
    integer(kind=shortInt), intent(in)  :: G
    self % G = G
    self % isMG = .true.
    khgjgj
  end subroutine makeMG

  subroutine buildCE(self,r,dir,E,w,Cell,Mat)
    class(particle), intent(out)                 :: self
    real(kind=defReal),dimension(3),intent(in)   :: r, dir
    real(kind=defReal),intent(in)                :: E, w
    integer(kind=shortInt),intent(in)            :: Cell, Mat

    self % r =r
    self % dir = dir
    self % E = E
    self % w = w
    self % CurrentCellIndex = Cell
    self % CurrentMaterialIndex = Mat

    self % isDead = .false.
    self % isMG = .false.
  end subroutine
    
  subroutine buildMG(self,r,dir,G,w,Cell,Mat)
    class(particle), intent(out)                 :: self
    real(kind=defReal),dimension(3),intent(in)   :: r, dir
    real(kind=defReal),intent(in)                :: w
    integer(kind=shortInt),intent(in)            :: G
    integer(kind=shortInt),intent(in)            :: Cell, Mat

    self % r =r
    self % dir = dir
    self % G = G
    self % w = w
    self % CurrentCellIndex = Cell
    self % CurrentMaterialIndex = Mat

    self % isDead = .false.
    self % isMG = .true.
  end subroutine

end module particle_class
