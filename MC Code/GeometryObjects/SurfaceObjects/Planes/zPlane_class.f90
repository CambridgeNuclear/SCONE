module zPlane_class

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
  character(nameLen),parameter :: TYPE_NAME    = 'zPlane'

  !!
  !! Constructor
  !!
  interface zPlane
    module procedure zPlane_fromDict
  end interface

  !!
  !! Simple plane perpendicular to z-axis
  !!
  type, public, extends(surface) :: zPlane
    real(defReal) :: z0 ! z-axis offset

  contains
    ! Initialisation & Indentification procedures
    procedure :: init
    procedure :: type
    procedure :: getDef

    ! Runtime procedures
    procedure :: evaluate
    procedure :: distance
    procedure :: normalVector
    procedure :: boundaryTransform

    procedure :: reflectAny

  end type zPlane

contains

  !!
  !! Initialise Z-Plane from offest
  !!
  subroutine init(self, z0, id)
    class(zPlane), intent(inout)  :: self
    real(defReal), intent(in)     :: z0
    integer(shortInt), intent(in) :: id

    self % z0 = z0
    call self % setId(id)

    ! Hash and store surface definition
    call self % hashSurfDef()

  end subroutine init

  !!
  !! Returns an initialised instance of zPlane for dictionary and name
  !!
  function zPlane_fromDict(dict) result(new)
    class(dictionary), intent(in)  :: dict
    type(zPlane)                   :: new
    integer(shortInt)              :: id
    real(defReal)                  :: z0
    character(100),parameter :: Here ='zPlane_fromDict ( zPlane_class.f90)'

    id = dict % getInt('id')
    if(id < 1) call fatalError(Here,'Invalid surface id provided')

    z0 = dict % getReal('z')
    call new % init(z0, id)

  end function zPlane_fromDict

  !!
  !! Evaluate plane distance
  !!
  elemental subroutine evaluate(self, res, r, shelf)
    class(zPlane), intent(in)      :: self
    real(defReal), intent(out)     :: res
    type(vector), intent(in)       :: r
    type(surfaceShelf), intent(in) :: shelf

    res = r % v(3) - self % z0

  end subroutine evaluate

  !!
  !! Return parameter character containing TYPE NAME
  !!
  elemental function type(self)
    class(zPlane), intent(in) :: self
    character(nameLen)        :: type

    type = TYPE_NAME

  end function type

  !!
  !! Return string with definition of this surface
  !!
  pure subroutine getDef(self,string)
    class(zPlane), intent(in)               :: self
    character(:),allocatable,intent(inout) :: string

    string = printSurfDef(TYPE_NAME, [self % z0])

  end subroutine getDef

  !!
  !! Calculate distance to plane along direction u
  !!
  elemental subroutine distance(self, dist, idx, r, u, shelf)
    class(zPlane), intent(in)        :: self
    real(defReal), intent(out)       :: dist
    integer(shortInt), intent(out)   :: idx
    type(vector), intent(in)         :: r
    type(vector),intent(in)          :: u
    type(surfaceShelf), intent(in)   :: shelf
    real(defReal)                    :: u3

    ! Set index
    idx = self % myIdx()

    u3   = u % v(3)
    dist = self % z0 - r % v(3)

    if ((u3 == ZERO) .OR. (abs(dist) < surface_tol)) then
      dist = INFINITY
      return
    end if

    dist = dist/u3
    if (dist < ZERO) dist = INFINITY
    return

  end subroutine distance

  !!
  !! Perform reflection
  !!
  elemental subroutine reflectAny(self, r, u, shelf)
    class(zPlane), intent(in)     :: self
    type(vector), intent(inout)   :: r
    type(vector), intent(inout)   :: u
    type(surfaceShelf),intent(in) :: shelf
    real(defReal)                 :: zDisplacement

    ! Reflect the particle z-coordinate across the plane
    zDisplacement = r % v(3) - self % z0
    r % v(3) = r % v(3) - TWO * zDisplacement

    ! Reflect the particle direction (independent of intersection point for plane)
    u % v(3) = -u % v(3)

  end subroutine reflectAny

  !!
  !! Returns vector normal to the plane
  !!
  elemental function normalVector(self, r, shelf) result(normal)
    class(zPlane), intent(in)      :: self
    type(vector), intent(in)       :: r
    type(surfaceShelf), intent(in) :: shelf
    type(vector)                   :: normal

    normal = [ZERO, ZERO, ONE]

  end function normalVector

  !!
  !! Apply boundary transformations
  !!
  subroutine boundaryTransform(self, r, u, isVacuum, shelf)
    class(zPlane), intent(in)      :: self
    type(vector), intent(inout)    :: r
    type(vector), intent(inout)    :: u
    logical(defBool), intent(out)  :: isVacuum
    type(surfaceShelf), intent(in) :: shelf
    character(100),parameter       :: Here ='boundaryTransform ( zPlane_class.f90)'

    call fatalError(Here,'zPlane does not support BC ')

  end subroutine boundaryTransform

end module zPlane_class
