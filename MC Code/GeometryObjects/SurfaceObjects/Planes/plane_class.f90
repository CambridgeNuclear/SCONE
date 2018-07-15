module plane_class
  use numPrecision
  use universalVariables
  use genericProcedures,  only : fatalError, dotProduct
  use vector_class,       only : vector
  use dictionary_class,   only : dictionary
  use surface_inter,      only : surface, printSurfDef, surfaceShelf

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

  end type plane

contains

  !!
  !! Initialise general plane from coefficients
  !!
  subroutine init(self, coeff, id)
    class(plane), intent(inout)             :: self
    real(defReal), dimension(4), intent(in) :: coeff
    integer(shortInt), intent(in)           :: id

    self % coeff = coeff
    call self % setId(id)

    ! Hash and store surface definition
    call self % hashSurfDef()

  end subroutine init

  !!
  !! Returns an initialised instance of plane for dictionary and name
  !!
  function plane_fromDict(dict) result(new)
    class(dictionary), intent(in)  :: dict
    type(plane)                    :: new
    integer(shortInt)              :: id
    real(defReal),dimension(4)     :: coeff
    character(100),parameter :: Here =' plane_fromDict ( plane_class.f90)'

    id = dict % getInt('id')
    if(id < 1) call fatalError(Here,'Invalid surface id provided')

    coeff = dict % getRealArray('coeff')
    call new % init(coeff, id)

  end function plane_fromDict

  !!
  !! Evaluate plane distance
  !!
  elemental subroutine evaluate(self,res, r)
    class(plane), intent(in)                :: self
    real(defReal), intent(out)              :: res
    type(vector), intent(in)                :: r

    res = r%v(1)*self%coeff(1) + r%v(2)*self%coeff(2) + r%v(3)*self%coeff(3) - self%coeff(4)

  end subroutine evaluate

  !!
  !! Return parameter character containing TYPE NAME
  !!
  elemental function type(self)
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
  elemental subroutine distance(self, dist, idx, r, u)
    class(plane), intent(in)          :: self
    real(defReal), intent(out)        :: dist
    integer(shortInt), intent(out)    :: idx
    type(vector), intent(in)          :: r
    type(vector), intent(in)          :: u
    real(defReal)                     :: denominator, numerator

    ! Set index
    idx = self % myIdx()

    denominator = self%coeff(1)*u%v(1) + self%coeff(2)*u%v(2) + self%coeff(3)*u%v(3)
    numerator = self%coeff(4) - self%coeff(1)*r%v(1) - self%coeff(2)*r%v(2) - self%coeff(3)*r%v(3)

    if ((denominator==ZERO) .OR. (abs(numerator) < surface_tol))  then
      dist = INFINITY
      return
    end if

    dist = numerator/denominator
    if (dist < ZERO) dist = INFINITY
    return

  end subroutine distance

  !!
  !! Perform reflection by the plane
  !!
  elemental subroutine reflect(self, r, u)
    class(plane), intent(in)       :: self
    type(vector), intent(inout)    :: r
    type(vector), intent(inout)    :: u
    type(vector)                   :: normal
    real(defReal)                  :: perpDist

    ! Translate the particle position across the plane
    call self % evaluate(perpDist, r)
    normal  = self % normalVector(r)
    r       = r - TWO * perpDist * normal

    ! Reflect the particle direction (independent of intersection point for plane)(assume normalised)
    u = u - TWO * (normal .dot. u) * normal

  end subroutine reflect

  !!
  !! Returns vector normal to the plane
  !!
  elemental function normalVector(self, r) result(normal)
    class(plane), intent(in)                :: self
    type(vector), intent(in)                :: r
    type(vector)                            :: normal

    normal = [self%coeff(1), self%coeff(2), self%coeff(3)]

  end function normalVector

end module plane_class
