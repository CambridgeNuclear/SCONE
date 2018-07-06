module plane_class
  use numPrecision
  use universalVariables
  use genericProcedures,  only : fatalError, dotProduct
  use dictionary_class,   only : dictionary
  use surface_inter,      only : surface, printSurfDef

  implicit none
  private

  !!
  !! Constants describing surface properties
  !!
  character(nameLen),parameter :: TYPE_NAME    = 'plane'

  !!
  !! Constructor
  !!
  interface plane
    module procedure plane_fromDict
  end interface

  !!
  !! General plane
  !!
  type, public, extends(surface) :: plane
    real(defReal), dimension(4), private :: coeff ! coefficients determining general plane (a,b,c,d)

  contains
    procedure :: init
    procedure :: evaluate
    procedure :: type
    procedure :: getDef
    procedure :: distanceToSurface
    procedure :: normalVector
    procedure :: whichSurface
    procedure :: boundaryTransform

    procedure,private :: reflectiveTransform

  end type plane

contains

  !!
  !! Initialise general plane from coefficients
  !!
  subroutine init(self, coeff, id, name)
    class(plane), intent(inout)             :: self
    real(defReal), dimension(4), intent(in) :: coeff
    integer(shortInt), intent(in), optional :: id
    character(*), optional, intent(in)      :: name

    self % coeff = coeff
    if(present(id)) self % id = id
    if(present(name)) self % name = name

  end subroutine init

  !!
  !! Returns an initialised instance of plane for dictionary and name
  !!
  function plane_fromDict(dict,name) result(new)
    class(dictionary), intent(in)  :: dict
    character(nameLen), intent(in) :: name
    type(plane)                    :: new
    integer(shortInt)              :: id
    real(defReal),dimension(4)     :: coeff
    character(100),parameter :: Here =' plane_fromDict ( plane_class.f90)'

    id = dict % getInt('id')
    if(id < 1) call fatalError(Here,'Invalid surface id provided')

    coeff = dict % getRealArray('coeff')
    call new % init(coeff, id, name)

  end function plane_fromDict

  !!
  !! Evaluate plane distance
  !!
  function evaluate(self, r) result(res)
    class(plane), intent(in)                :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: res

    res = r(1)*self%coeff(1) + r(2)*self%coeff(2) + r(3)*self%coeff(3) - self%coeff(4)

  end function evaluate

  !!
  !! Return parameter character containing TYPE NAME
  !!
  function type(self)
    class(plane), intent(in) :: self
    character(nameLen)     :: type

    type = TYPE_NAME

  end function type

  !!
  !! Return string with definition of this surface
  !!
  pure subroutine getDef(self,string)
    class(plane), intent(in)               :: self
    character(:),allocatable,intent(inout) :: string

    string = printSurfDef(TYPE_NAME, self % coeff)

  end subroutine getDef

  !!
  !! Calculate distance to plane along direction u
  !!
  function distanceToSurface(self,r,u) result(distance)
    class(plane), intent(in)                :: self
    real(defReal), dimension(3), intent(in) :: r, u
    real(defReal)                           :: distance, denominator, numerator

    denominator = self%coeff(1)*u(1) + self%coeff(2)*u(2) + self%coeff(3)*u(3)
    numerator = self%coeff(4) - self%coeff(1)*r(1) - self%coeff(2)*r(2) - self%coeff(3)*r(3)

    if ((denominator==ZERO) .OR. (abs(numerator) < surface_tol))  then
      distance = INFINITY
      return
    end if

    distance = numerator/denominator
    if (distance < ZERO) distance = INFINITY
    return

  end function distanceToSurface

  !!
  !! Perform reflection
  !!
  subroutine reflectiveTransform(self,r,u)
    class(plane), intent(in)                   :: self
    real(defReal), dimension(3), intent(inout) :: r, &
                                                  u
    real(defReal), dimension(3)                :: normal
    real(defReal)                              :: magSquared, perpDistance

    ! Re-examine whether this should be commented after defining the transport operator
    !p%r = p%r + p%dir*distance

    ! Translate the particle position across the plane
    normal = self%normalVector(r)
    magSquared = dotProduct(normal,normal)
    perpDistance = self%evaluate(r)
    r = r - TWO*perpDistance*normal

    ! Reflect the particle direction (independent of intersection point for plane)(assume normalised)
    u = u - TWO*dotProduct(normal,u)*normal

  end subroutine reflectiveTransform

  !!
  !! Returns vector normal to the plane
  !!
  function normalVector(self,r) result(normal)
    class(plane), intent(in)                :: self
    real(defReal), dimension(3)             :: normal
    real(defReal), dimension(3), intent(in) :: r

    normal = [self%coeff(1), self%coeff(2), self%coeff(3)]

  end function normalVector

  !!
  !! Give an error: this routine should not be called for a non-compound surface
  !!
  function whichSurface(self, r, u) result(surfPointer)
    class(plane), intent(in)                :: self
    real(defReal), dimension(3), intent(in) :: r, u
    class(surface), pointer                 :: surfPointer
    character(100),parameter :: Here = 'whichSurface (plane_class.f90)'

    call fatalError(Here,'This function should never be called for a simple surface')

  end function whichSurface

  !!
  !! Apply boundary transformation
  !!
  subroutine boundaryTransform(self, r, u, isVacuum)
    class(plane), intent(in)                   :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u
    logical(defBool), intent(inout)            :: isVacuum
    character(100),parameter :: Here = 'boundaryTransform (plane_class.f90)'

    call fatalError(Here,'This surface does not support boundary conditions')

  end subroutine boundaryTransform

end module plane_class
