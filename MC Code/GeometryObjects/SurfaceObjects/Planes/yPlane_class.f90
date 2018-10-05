module yPlane_class

  use numPrecision
  use universalVariables
  use genericProcedures,  only : fatalError
  use vector_class,       only : vector
  use dictionary_class,   only : dictionary
  use surface_inter,      only : surface, printSurfDef, surfaceShelf

  implicit none
  private

  !!
  !! Constants describing surface properties
  !!
  character(nameLen),parameter :: TYPE_NAME    = 'yPlane'

  !!
  !! Constructor
  !!
  interface yPlane
    module procedure yPlane_fromDict
  end interface

  !!
  !! Simple plane perpendicular to y-axis
  !!
  type, public, extends(surface) :: yPlane
    real(defReal) :: y0 ! y-axis offset

  contains
    ! Initialisation & Indentification procedures
    procedure :: init
    procedure :: type
    procedure :: getDef
    procedure :: boundingBox

    ! Runtime procedures
    procedure :: evaluate
    procedure :: distance
    procedure :: normalVector

    ! Type specific procedures
    procedure :: reflect

  end type yPlane


contains

  !!
  !! Initialise Y-Plane from offset
  !!
  subroutine init(self, y0, id)
    class(yPlane), intent(inout)            :: self
    real(defReal), intent(in)               :: y0
    integer(shortInt), intent(in), optional :: id

    self % y0 = y0
    call self % setId(id)

    ! Hash and store surface definition
    call self % hashSurfDef()

  end subroutine init

  !!
  !! Returns an initialised instance of yPlane for dictionary and name
  !!
  function yPlane_fromDict(dict) result(new)
    class(dictionary), intent(in)  :: dict
    type(yPlane)                   :: new
    integer(shortInt)              :: id
    real(defReal)                  :: y0
    character(100),parameter :: Here ='yPlane_fromDict ( yPlane_class.f90)'

    call dict % get(id, 'id')
    if(id < 1) call fatalError(Here,'Invalid surface id provided')

    call dict % get(y0, 'y')
    call new % init(y0, id)


  end function yPlane_fromDict

  !!
  !! Evaluate plane distance
  !!
  elemental subroutine evaluate(self, res, r)
    class(yPlane), intent(in)     :: self
    real(defReal), intent(out)    :: res
    type(vector), intent(in)      :: r

    res = r % v(2) - self % y0

  end subroutine evaluate

  !!
  !! Return parameter character containing TYPE NAME
  !!
  elemental function type(self)
    class(yPlane), intent(in) :: self
    character(nameLen)        :: type

    type = TYPE_NAME

  end function type

  !!
  !! Return string with definition of this surface
  !!
  pure subroutine getDef(self,string)
    class(yPlane), intent(in)              :: self
    character(:),allocatable,intent(inout) :: string

    string = printSurfDef(TYPE_NAME, [self % y0])

  end subroutine getDef

  !!
  !! Returns an axis alligned bouding box of surface -ve halfspace
  !!
  pure subroutine boundingBox(self,origin, halfwidth)
    class(yPlane), intent(in)               :: self
    real(defReal), dimension(3),intent(out) :: origin
    real(defReal), dimension(3),intent(out) :: halfwidth

    origin    = ZERO
    halfwidth = INFINITY

    origin(2)      = (INFINITY + self % y0) * HALF
    halfwidth(2)   = (INFINITY + self % y0) * HALF

  end subroutine boundingBox

  !!
  !! Calculate distance to plane along direction u
  !!
  elemental subroutine distance(self, dist, idx, r, u)
    class(yPlane), intent(in)       :: self
    real(defReal), intent(out)      :: dist
    integer(shortInt), intent(out)  :: idx
    type(vector), intent(in)        :: r
    type(vector), intent(in)        :: u
    real(defReal)                   :: u2

    ! Set index
    idx = self % myIdx()

    u2   = u % v(2)
    dist = self % y0 - r % v(2)

    if ((u2 == ZERO) .OR. (abs(dist) < surface_tol)) then
      dist = INFINITY
      return
    end if

    dist = dist/u2
    if (dist < ZERO) dist = INFINITY
    return

  end subroutine distance

  !!
  !! Perform reflection by the plane
  !!
  elemental subroutine reflect(self, r, u)
    class(yPlane), intent(in)        :: self
    type(vector), intent(inout)      :: r
    type(vector), intent(inout)      :: u
    real(defReal)                    :: yDisplacement

    ! Reflect the particle y-coordinate across the plane
    yDisplacement = r % v(2) - self % y0
    r % v(2) = r % v(2) - TWO * yDisplacement

    ! Reflect the particle direction (independent of intersection point for plane)
    u % v(2) = -u % v(2)

  end subroutine reflect

  !!
  !! Returns vector normal to the plane
  !!
  elemental function normalVector(self, r) result(normal)
    class(yPlane), intent(in)      :: self
    type(vector), intent(in)       :: r
    type(vector)                   :: normal

    normal = [ZERO, ONE, ZERO]

  end function normalVector

end module yPlane_class
