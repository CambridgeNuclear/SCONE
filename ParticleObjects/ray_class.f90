module ray_class

  use numPrecision
  use universalVariables
  use genericProcedures
  use coord_class,       only : coordList
  use RNG_class,         only : RNG

  implicit none
  private

  !!
  !! Ray to be used with Random Ray Method implementations
  !!
  !! Very simple version of a particle, simply used for traversing geometry
  !!
  !! Public Members:
  !!   coords   -> Position information of the ray
  !!
  type, public :: ray
      
    ! Ray phase space data
    type(coordList)                          :: coords  ! Space/angle co-ordinates

    ! Ray processing information
    class(RNG), pointer        :: pRNG  => null()  ! Pointer to RNG associated with the ray
    integer(shortInt)          :: geomIdx          ! Index of the geometry used by the ray 

  contains
     ! Build procedures
    procedure                  :: build 
    procedure                  :: kill

    ! Inquiry about coordinates
    procedure                  :: rLocal
    procedure                  :: rGlobal
    procedure                  :: dirLocal
    procedure                  :: dirGlobal
    procedure                  :: nesting
    procedure                  :: getCellIdx
    procedure                  :: getUniIdx
    procedure                  :: matIdx

    ! Operations on coordinates
    procedure            :: moveGlobal
    procedure            :: moveLocal
    procedure            :: rotate
    procedure            :: teleport
    procedure            :: point
    procedure            :: takeAboveGeom
    procedure            :: setMatIdx
    
    ! Debug procedures
    procedure            :: display

  end type ray

contains

!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Particle build and assignment procedures
!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Initialise ray
  !!   r   -> Global Position
  !!   dir -> Global direction
  !!   ng  -> Number of energy groups
  !!
  subroutine build(self, r, dir)
    class(ray), intent(inout)               :: self
    real(defReal),dimension(3),intent(in)   :: r
    real(defReal),dimension(3),intent(in)   :: dir
    character(nameLen), parameter           :: Here = 'build (ray_class.f90)'

    call self % coords % init(r, dir)

  end subroutine build

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(ray), intent(inout) :: self

    call self % coords % kill()

  end subroutine kill

!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Particle coordinates inquiry procedures
!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Return the position either at the deepest nested level or a specified level
  !!
  function rLocal(self,n)result(r)
    class(ray), intent(in)                  :: self
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
    class(ray), intent(in)      :: self
    real(defReal), dimension(3) :: r

    r = self % coords % lvl(1) % r

  end function rGlobal

  !!
  !! Return the direction either at the deepest nested level or at a specified level
  !!
  function dirLocal(self,n)result(dir)
    class(ray), intent(in)      :: self
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
    class(ray), intent(in)      :: self
    real(defReal), dimension(3) :: dir

    dir = self % coords % lvl(1) % dir

  end function dirGlobal

  !!
  !! Return the lowest nesting level of the ray
  !!
  function nesting(self) result(n)
    class(ray), intent(in)      :: self
    integer(shortInt)           :: n

    n = self % coords % nesting

  end function nesting

  !!
  !! Return cell index at a given nesting level n
  !!  If no n is given return for lowest nesting level
  !!
  function getCellIdx(self,n) result(idx)
    class(ray), intent(in)                  :: self
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
  function getUniIdx(self,n) result(idx)
    class(ray), intent(in)                  :: self
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
    class(ray), intent(in) :: self
    integer(shortInt)      :: idx

    idx = self % coords % matIdx

  end function matIdx

!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Ray operations on coordinates procedures
!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Move the ray above the geometry
  !! NOTE: regionID & matIdx will be reset!!!
  !!
  subroutine moveGlobal(self,distance)
    class(ray), intent(inout) :: self
    real(defReal), intent(in) :: distance

    call self % coords % moveGlobal(distance)

  end subroutine moveGlobal

  !!
  !! Move ray in local co-ordinates down to nesting level n
  !!
  subroutine moveLocal(self,distance,n)
    class(ray), intent(inout)      :: self
    real(defReal), intent(in)      :: distance
    integer(shortInt), intent(in)  :: n

    call self % coords % moveLocal(distance,n)

  end subroutine moveLocal

  !!
  !! Rotate ray
  !!  mu  -> cosine of deflection from current direction
  !!  phi -> azimuthal angle of rotation
  !!
  subroutine rotate(self,mu,phi)
    class(ray), intent(inout) :: self
    real(defReal), intent(in) :: mu
    real(defReal), intent(in) :: phi

    call self % coords % rotate(mu,phi)

  end subroutine rotate

  !!
  !! Place ray at an arbitrary point in above the geometry
  !!
  subroutine teleport(self, r)
    class(ray), intent(inout)               :: self
    real(defReal), dimension(3), intent(in) :: r

    call self % coords % assignPosition(r)

  end subroutine teleport

  !!
  !! Point ray in direction dir in highest nesting level
  !! Propagates new direction to lower levels
  !!
  subroutine point(self, dir)
    class(ray), intent(inout)               :: self
    real(defReal), dimension(3), intent(in) :: dir

    call self % coords % assignDirection(dir)

  end subroutine point

  !!
  !! Resets the ray's nesting level
  !!
  subroutine takeAboveGeom(self)
    class(ray), intent(inout) :: self

    call self % coords % takeAboveGeom()

  end subroutine takeAboveGeom

  !!
  !! Set Material index for testing purposes
  !!
  pure subroutine setMatIdx(self,matIdx)
    class(ray), intent(inout)      :: self
    integer(shortInt), intent(in)  :: matIdx

    self % coords % matIdx = matIdx

  end subroutine setMatIdx

!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Ray debug procedures
!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Display state of a ray
  !!
  subroutine display(self)
    class(ray), intent(in) :: self

    print *, self % coords % matIdx

  end subroutine display


end module ray_class
