!
! Module containing the 3 planes perpendicular to the x, y, or z axis and the general plane
!
module plane_class
  use numPrecision
  use universalVariables
  use genericProcedures,    only : fatalError, dotProduct

  use surface_inter,        only : surface

  implicit none
  private

  !!
  !! General plane
  !!
  type, public, extends(surface) :: plane
    real(defReal), dimension(4), private :: coeff ! coefficients determining general plane (a,b,c,d)

  contains
    procedure :: init
    procedure :: evaluate
    procedure :: distanceToSurface
    procedure :: reflectiveTransform
    procedure :: normalVector
    procedure :: whichSurface
    procedure :: setBoundaryConditions
    procedure :: boundaryTransform

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
  !! Evaluate plane distance
  !!
  function evaluate(self, r) result(res)
    class(plane), intent(in)                :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: res

    res = r(1)*self%coeff(1) + r(2)*self%coeff(2) + r(3)*self%coeff(3) - self%coeff(4)

  end function evaluate

  !!
  !! Calculate distance to plane along direction u
  !!
  function distanceToSurface(self,r,u) result(distance)
    class(plane), intent(in)                :: self
    real(defReal), dimension(3), intent(in) :: r, &
                                               u
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
    return

  end function normalVector

  !!
  !! Give an error: this routine should not be called for a non-compound surface
  !!
  function whichSurface(self, r, u) result(surfPointer)
    class(plane), intent(in)                :: self
    real(defReal), dimension(3), intent(in) :: r, u
    class(surface), pointer                 :: surfPointer

    call fatalError('whichPlane','This function should never be called for a simple surface')

  end function whichSurface

  !!
  !! Crash on attempting to set boundary conditions for a plane object
  !!
  subroutine setBoundaryConditions(self, BC)
    class(plane), intent(inout)                 :: self
    integer(shortInt), dimension(6), intent(in) :: BC
    call fatalError('setBoundaryConditionsPlane','Boundary conditions may not be set for a plane surface')
  end subroutine setBoundaryConditions

  !!
  !! Apply boundary transformation
  !!
  subroutine boundaryTransform(self, r, u, isVacuum)
    class(plane), intent(in)                   :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u
    logical(defBool), intent(inout)            :: isVacuum

    call fatalError('boundaryTransform, plane','This surface does not support boundary conditions')

  end subroutine boundaryTransform

end module plane_class
