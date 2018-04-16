!
! Module containing the 3 planes perpendicular to the x, y, or z axis and the general plane
!
module plane_class
  use numPrecision
  use genericProcedures
  use universalVariables

  use surface_class

  implicit none
  private

  type, public, extends(surface) :: xPlane
    real(defReal) :: x0 ! x-axis offset
  contains
    procedure :: init => initXPlane
    procedure :: evaluate => evaluateXPlane
    procedure :: distanceToSurface => distanceToXPlane
    procedure :: reflectiveTransform => reflectiveTransformXPlane
    procedure :: normalVector => normalVectorXPlane
    procedure :: whichSurface => whichXPlane
    procedure :: setBoundaryConditions => setBoundaryConditionsXPlane
  end type xPlane

  type, public, extends(surface) :: yPlane
    real(defReal) :: y0 ! x-axis offset
  contains
    procedure :: init => initYPlane
    procedure :: evaluate => evaluateYPlane
    procedure :: distanceToSurface => distanceToYPlane
    procedure :: reflectiveTransform => reflectiveTransformYPlane
    procedure :: normalVector => normalVectorYPlane
    procedure :: whichSurface => whichYPlane
    procedure :: setBoundaryConditions => setBoundaryConditionsYPlane
  end type yPlane

  type, public, extends(surface) :: zPlane
    real(defReal) :: z0 ! x-axis offset
  contains
    procedure :: init => initZPlane
    procedure :: evaluate => evaluateZPlane
    procedure :: distanceToSurface => distanceToZPlane
    procedure :: reflectiveTransform => reflectiveTransformZPlane
    procedure :: normalVector => normalVectorZPlane
    procedure :: whichSurface => whichZPlane
    procedure :: setBoundaryConditions => setBoundaryConditionsZPlane
  end type zPlane

  type, public, extends(surface) :: plane
    real(defReal), dimension(4), private :: coeff ! coefficients determining general plane
  contains
    procedure :: init => initPlane
    procedure :: evaluate => evaluatePlane
    procedure :: distanceToSurface => distanceToPlane
    procedure :: reflectiveTransform => reflectiveTransformPlane
    procedure :: normalVector => normalVectorPlane
    procedure :: whichSurface => whichPlane
    procedure :: setBoundaryConditions => setBoundaryConditionsPlane
  end type plane

contains

!
! X-Plane procedures
!
  subroutine initXPlane(self, x0, id, name)
    class(xPlane), intent(inout) :: self
    real(defReal), intent(in) :: x0
    integer(shortInt), intent(in), optional :: id
    character(*), optional, intent(in) :: name

    self % x0 = x0
    if(present(id)) self % id = id
    if(present(name)) self % name = name

  end subroutine initXPlane

  function evaluateXPlane(self, r) result(res)
    implicit none
    class(xPlane), intent(in) :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal) :: res

    res = r(1) - self % x0

  end function evaluateXPlane

  function distanceToXPlane(self,r,u)result(distance)
    implicit none
    class(xPlane), intent(in) :: self
    real(defReal), dimension(3), intent(in) :: r, &
                                               u
    real(defReal) :: distance, u1

    u1 = u(1)
    distance = self%x0 - r(1)
    if ((u1==0.0) .OR. (abs(distance) < surface_tol)) then
      distance = INFINITY
      return
    end if

    distance = distance/u1
    if (distance < 0.0) distance = INFINITY
    return

  end function distanceToXPlane

  subroutine reflectiveTransformXPlane(self, r, u)
    implicit none
    class(xPlane), intent(in) :: self
    real(defReal), dimension(3), intent(inout) :: r, &
                                                  u
    real(defReal) :: xDisplacement

    ! Reflect the particle x-coordinate across the plane
    xDisplacement = r(1) - self%x0
    r(1) = r(1) - 2*xDisplacement

    ! Reflect the particle direction (independent of intersection point for plane)
    u(1) = -u(1)

  end subroutine reflectiveTransformXPlane

  function normalVectorXPlane(self,r)result(normal)
    implicit none
    class(xPlane), intent(in) :: self
    real(defReal), dimension(3) :: normal
    real(defReal), dimension(3), intent(in) :: r

    normal = [1, 0, 0]
    return

  end function normalVectorXPlane

  !!
  !! Give an error: this routine should not be called for a non-compound surface
  !!
  function whichXPlane(self, r, u) result(surfPointer)
    class(xPlane), intent(in) :: self
    real(defReal), dimension(3), intent(in) :: r, u
    class(surface), pointer :: surfPointer
    call fatalError('whichXPlane','This function should never be called for a simple surface')
  end function whichXPlane

  !!
  !! Crash on attempting to set boundary conditions for a plane object
  !!
  subroutine setBoundaryConditionsXPlane(self, BC)
    class(xPlane), intent(inout) :: self
    integer(shortInt), dimension(6), intent(in) :: BC
    call fatalError('setBoundaryConditionsXPlane','Boundary conditions may not be set for a plane surface')
  end subroutine setBoundaryConditionsXPlane

!
! Y-Plane procedures
!
  subroutine initYPlane(self, y0, id, name)
    implicit none
    class(yPlane), intent(inout) :: self
    real(defReal), intent(in) :: y0
    integer(shortInt), intent(in), optional :: id
    character(*), optional, intent(in) :: name

    self % y0 = y0
    if(present(id)) self % id = id
    if(present(name)) self % name = name

  end subroutine initYPlane

  function evaluateYPlane(self, r) result(res)
    implicit none
    class(yPlane), intent(in) :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal) :: res

    res = r(2) - self % y0

  end function evaluateYPlane

  function distanceToYPlane(self,r,u)result(distance)
    implicit none
    class(yPlane), intent(in) :: self
    real(defReal), dimension(3), intent(in) :: r, &
                                               u
    real(defReal) :: distance, v

    v = u(2)
    distance = self%y0 - r(2)
    if ((v==0.0) .OR. (abs(distance) < surface_tol)) then
      distance = INFINITY
      return
    end if

    distance = distance/v
    if (distance < 0.0) distance = INFINITY
    return

  end function distanceToYPlane

  subroutine reflectiveTransformYPlane(self,r,u)
    implicit none
    class(yPlane), intent(in) :: self
    real(defReal), dimension(3), intent(inout) :: r, &
                                                  u
    real(defReal) :: yDisplacement

    ! Reflect the particle y-coordinate across the plane
    yDisplacement = r(2) - self%y0
    r(2) = r(2) - 2*yDisplacement

    ! Reflect the particle direction (independent of intersection point for plane)
    u(2) = -u(2)

  end subroutine reflectiveTransformYPlane

  function normalVectorYPlane(self,r)result(normal)
    implicit none
    class(yPlane), intent(in) :: self
    real(defReal), dimension(3) :: normal
    real(defReal), dimension(3), intent(in) :: r

    normal = [0, 1, 0]
    return

  end function normalVectorYPlane

  !
  ! Give an error: this routine should not be called for a non-compound surface
  !
  function whichYPlane(self, r, u) result(surfPointer)
    implicit none
    class(yPlane), intent(in) :: self
    real(defReal), dimension(3), intent(in) :: r, u
    class(surface), pointer :: surfPointer
    call fatalError('whichYPlane','This function should never be called for a simple surface')
  end function whichYPlane

  !!
  !! Crash on attempting to set boundary conditions for a plane object
  !!
  subroutine setBoundaryConditionsYPlane(self, BC)
    class(yPlane), intent(inout) :: self
    integer(shortInt), dimension(6), intent(in) :: BC
    call fatalError('setBoundaryConditionsYPlane','Boundary conditions may not be set for a plane surface')
  end subroutine setBoundaryConditionsYPlane

!
! Z-Plane procedures
!
  subroutine initZPlane(self, z0, id, name)
    implicit none
    class(zPlane), intent(inout) :: self
    real(defReal), intent(in) :: z0
    integer(shortInt), intent(in), optional :: id
    character(*), optional, intent(in) :: name

    self % z0 = z0
    if(present(id)) self % id = id
    if(present(name)) self % name = name

  end subroutine initZPlane

  function evaluateZPlane(self, r) result(res)
    implicit none
    class(zPlane), intent(in) :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal) :: res

    res = r(3) - self % z0

  end function evaluateZPlane

  function distanceToZPlane(self,r,u)result(distance)
    implicit none
    class(zPlane), intent(in) :: self
    real(defReal), dimension(3), intent(in) :: r, &
                                               u
    real(defReal) :: distance, w

    w = u(3)
    distance = self%z0 - r(3)
    if ((w==0.0) .OR. (abs(distance) < surface_tol)) then
      distance = INFINITY
      return
    end if

    distance = (self % z0 - r(3))/w
    if (distance < 0.0) distance = INFINITY
    return

  end function distanceToZPlane

  subroutine reflectiveTransformZPlane(self,r,u)
    implicit none
    class(zPlane), intent(in) :: self
    real(defReal), dimension(3), intent(inout) :: r, &
                                                  u
    real(defReal) :: zDisplacement

    ! Reflect the particle z-coordinate across the plane
    zDisplacement = r(3) - self%z0
    r(3) = r(3) - 2*zDisplacement

    ! Reflect the particle direction (independent of intersection point for plane)
    u(3) = -u(3)

  end subroutine reflectiveTransformZPlane

  function normalVectorZPlane(self,r)result(normal)
    implicit none
    class(zPlane), intent(in) :: self
    real(defReal), dimension(3) :: normal
    real(defReal), dimension(3), intent(in) :: r

    normal = [0, 0, 1]
    return

  end function normalVectorZPlane

  !
  ! Give an error: this routine should not be called for a non-compound surface
  !
  function whichZPlane(self, r, u) result(surfPointer)
    implicit none
    class(zPlane), intent(in) :: self
    real(defReal), dimension(3), intent(in) :: r, u
    class(surface), pointer :: surfPointer
    call fatalError('whichZPlane','This function should never be called for a simple surface')
  end function whichZPlane

  !!
  !! Crash on attempting to set boundary conditions for a plane object
  !!
  subroutine setBoundaryConditionsZPlane(self, BC)
    class(zPlane), intent(inout) :: self
    integer(shortInt), dimension(6), intent(in) :: BC
    call fatalError('setBoundaryConditionsZPlane','Boundary conditions may not be set for a plane surface')
  end subroutine setBoundaryConditionsZPlane

!
! General Plane procedures
!
  subroutine initPlane(self, coeff, id, name)
    implicit none
    class(plane), intent(inout) :: self
    real(defReal), dimension(4), intent(in) :: coeff
    integer(shortInt), intent(in), optional :: id
    character(*), optional, intent(in) :: name

    self % coeff = coeff
    if(present(id)) self % id = id
    if(present(name)) self % name = name

  end subroutine initPlane

  function evaluatePlane(self, r) result(res)
    implicit none
    class(plane), intent(in) :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal) :: res

    res = r(1)*self%coeff(1) + r(2)*self%coeff(2) + r(3)*self%coeff(3) - self%coeff(4)

  end function evaluatePlane

  function distanceToPlane(self,r,u)result(distance)
    implicit none
    class(plane), intent(in) :: self
    real(defReal), dimension(3), intent(in) :: r, &
                                               u
    real(defReal) :: distance, denominator, numerator

    denominator = self%coeff(1)*u(1) + self%coeff(2)*u(2) + self%coeff(3)*u(3)
    numerator = self%coeff(4) - self%coeff(1)*r(1) - self%coeff(2)*r(2) - self%coeff(3)*r(3)
    if ((denominator==0.0) .OR. (abs(numerator) < surface_tol))  then
      distance = INFINITY
      return
    end if

    distance = numerator/denominator
    if (distance < 0.0) distance = INFINITY
    return

  end function distanceToPlane

  subroutine reflectiveTransformPlane(self,r,u)
    implicit none
    class(plane), intent(in) :: self
    real(defReal), dimension(3), intent(inout) :: r, &
                                                  u
    real(defReal), dimension(3) :: normal
    real(defReal) :: magSquared, perpDistance

    ! Re-examine whether this should be commented after defining the transport operator
    !p%r = p%r + p%dir*distance

    ! Translate the particle position across the plane
    normal = self%normalVector(r)
    magSquared = dotProduct(normal,normal)
    perpDistance = self%evaluate(r)
    r = r - 2*perpDistance*normal

    ! Reflect the particle direction (independent of intersection point for plane)(assume normalised)
    u = u - 2*dotProduct(normal,u)*normal

  end subroutine reflectiveTransformPlane

  function normalVectorPlane(self,r)result(normal)
    implicit none
    class(plane), intent(in) :: self
    real(defReal), dimension(3) :: normal
    real(defReal), dimension(3), intent(in) :: r

    normal = [self%coeff(1), self%coeff(2), self%coeff(3)]
    return

  end function normalVectorPlane

  !
  ! Give an error: this routine should not be called for a non-compound surface
  !
  function whichPlane(self, r, u) result(surfPointer)
    implicit none
    class(plane), intent(in) :: self
    real(defReal), dimension(3), intent(in) :: r, u
    class(surface), pointer :: surfPointer
    call fatalError('whichPlane','This function should never be called for a simple surface')
  end function whichPlane

  !!
  !! Crash on attempting to set boundary conditions for a plane object
  !!
  subroutine setBoundaryConditionsPlane(self, BC)
    class(plane), intent(inout) :: self
    integer(shortInt), dimension(6), intent(in) :: BC
    call fatalError('setBoundaryConditionsPlane','Boundary conditions may not be set for a plane surface')
  end subroutine setBoundaryConditionsPlane

end module plane_class
