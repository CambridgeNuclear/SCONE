module xPlane_class

  use numPrecision
  use universalVariables
  use genericProcedures,    only : fatalError
  use vector_class,         only : vector
  use dictionary_class,     only : dictionary
  use surface_inter,        only : surface, printSurfDef, surfaceShelf

  implicit none
  private

  !!
  !! Constants describing surface properties
  !!
  character(nameLen),parameter :: TYPE_NAME    = 'xPlane'

  !!
  !! Constructor
  !!
  interface xPlane
    module procedure xPlane_fromDict
  end interface

  !!
  !! Simple plane perpendicular to x-axis
  !!
  type, public, extends(surface) :: xPlane
    real(defReal) :: x0 ! x-axis offset

  contains
    ! Initialisation & Indentification procedures
    procedure :: init
    procedure :: type
    procedure :: getDef

    ! Runtime procedures
    procedure :: evaluate
    procedure :: distance
    procedure :: normalVector

    ! Type specific procedures
    procedure :: reflect

  end type xPlane

contains

  !!
  !! Initialise X-Plane from offset
  !!
  subroutine init(self, x0, id)
    class(xPlane), intent(inout)  :: self
    real(defReal), intent(in)     :: x0
    integer(shortInt), intent(in) :: id

    self % x0 = x0
    call self % setId(id)

    ! Hash and store surface definition
    call self % hashSurfDef()

  end subroutine init

  !!
  !! Returns an initialised instance of xPlane for dictionary and name
  !!
  function xPlane_fromDict(dict) result(new)
    class(dictionary), intent(in)  :: dict
    type(xPlane)                   :: new
    integer(shortInt)              :: id
    real(defReal)                  :: x0
    character(100),parameter :: Here =' xPlane_fromDict ( xPlane_class.f90)'

    id = dict % getInt('id')
    if(id < 1) call fatalError(Here,'Invalid surface id provided')

    x0 = dict % getReal('x')
    call new % init(x0, id)

  end function xPlane_fromDict

  !!
  !!  Evaluate plane distance
  !!
  elemental subroutine evaluate(self, res, r)
    class(xPlane), intent(in)      :: self
    real(defReal), intent(out)     :: res
    type(vector), intent(in)       :: r

    res = r % v(1) - self % x0

  end subroutine evaluate

  !!
  !! Return parameter character containing TYPE NAME
  !!
  elemental function type(self)
    class(xPlane), intent(in) :: self
    character(nameLen)        :: type

    type = TYPE_NAME

  end function type

  !!
  !! Return string with definition of this surface
  !!
  pure subroutine getDef(self,string)
    class(xPlane), intent(in)              :: self
    character(:),allocatable,intent(inout) :: string

    string = printSurfDef(TYPE_NAME, [self % x0])

  end subroutine getDef

  !!
  !! Calculate distance to plane along direction u
  !!
  elemental subroutine distance(self, dist, idx, r, u)
    class(xPlane), intent(in)      :: self
    real(defReal), intent(out)     :: dist
    integer(shortInt), intent(out) :: idx
    type(vector), intent(in)       :: r
    type(vector), intent(in)       :: u
    real(defReal)                  :: u1

    ! Set index
    idx = self % myIdx()

    u1   = u % v(1)
    dist = self % x0 - r % v(1)

    if ((u1==ZERO) .OR. (abs(dist) < surface_tol)) then
      dist = INFINITY
      return
    end if

    dist = dist/u1
    if (dist < ZERO) dist = INFINITY

  end subroutine distance

  !!
  !! Perform reflection by the plane
  !!
  elemental subroutine reflect(self, r, u)
    class(xPlane), intent(in)      :: self
    type(vector), intent(inout)    :: r
    type(vector), intent(inout)    :: u
    real(defReal)                  :: xDisplacement

    ! Reflect the particle x-coordinate across the plane
    xDisplacement = r % v(1) - self % x0
    r % v(1) = r % v(1) - TWO * xDisplacement

    ! Reflect the particle direction (independent of intersection point for plane)
    u % v(1) = -u % v(1)

  end subroutine reflect

  !!
  !! Returns vector normal to the plane
  !!
  elemental function normalVector(self, r) result(normal)
    class(xPlane), intent(in)      :: self
    type(vector), intent(in)       :: r
    type(vector)                   :: normal

    normal = [ONE, ZERO, ZERO]

  end function normalVector
end module xPlane_class
