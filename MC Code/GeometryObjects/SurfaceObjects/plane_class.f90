!
! Module containing the 3 planes perpendicular to the x, y, or z axis and the general plane
!
module plane_class
  use numPrecision
  use genericProcedures, only : fatalError, dotProduct
  use universalVariables

  use surface_inter

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
    procedure :: boundaryTransform => boundaryTransformXPlane
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
    procedure :: boundaryTransform => boundaryTransformYPlane
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
    procedure :: boundaryTransform => boundaryTransformZPlane
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
    procedure :: boundaryTransform => boundaryTransformPlane
  end type plane

contains

!!
!! X-Plane procedures
!!
  subroutine initXPlane(self, x0, id, name)
    class(xPlane), intent(inout)            :: self
    real(defReal), intent(in)               :: x0
    integer(shortInt), intent(in), optional :: id
    character(*), optional, intent(in)      :: name

    self % x0 = x0
    if(present(id)) self % id = id
    if(present(name)) self % name = name

  end subroutine initXPlane

  function evaluateXPlane(self, r) result(res)
    class(xPlane), intent(in)               :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: res

    res = r(1) - self % x0

  end function evaluateXPlane

  function distanceToXPlane(self,r,u)result(distance)
    class(xPlane), intent(in)               :: self
    real(defReal), dimension(3), intent(in) :: r, &
                                               u
    real(defReal)                           :: distance, u1

    u1 = u(1)
    distance = self%x0 - r(1)
    if ((u1==ZERO) .OR. (abs(distance) < surface_tol)) then
      distance = INFINITY
      return
    end if

    distance = distance/u1
    if (distance < ZERO) distance = INFINITY
    return

  end function distanceToXPlane

  subroutine reflectiveTransformXPlane(self, r, u)
    class(xPlane), intent(in)                  :: self
    real(defReal), dimension(3), intent(inout) :: r, &
                                                  u
    real(defReal)                              :: xDisplacement

    ! Reflect the particle x-coordinate across the plane
    xDisplacement = r(1) - self%x0
    r(1) = r(1) - TWO*xDisplacement

    ! Reflect the particle direction (independent of intersection point for plane)
    u(1) = -u(1)

  end subroutine reflectiveTransformXPlane

  function normalVectorXPlane(self,r)result(normal)
    class(xPlane), intent(in)               :: self
    real(defReal), dimension(3)             :: normal
    real(defReal), dimension(3), intent(in) :: r

    normal = [ONE, ZERO, ZERO]
    return

  end function normalVectorXPlane

  !!
  !! Give an error: this routine should not be called for a non-compound surface
  !!
  function whichXPlane(self, r, u) result(surfPointer)
    class(xPlane), intent(in)               :: self
    real(defReal), dimension(3), intent(in) :: r, u
    class(surface), pointer                 :: surfPointer
    call fatalError('whichXPlane','This function should never be called for a simple surface')
  end function whichXPlane

  !!
  !! Crash on attempting to set boundary conditions for a plane object
  !!
  subroutine setBoundaryConditionsXPlane(self, BC)
    class(xPlane), intent(inout)                :: self
    integer(shortInt), dimension(6), intent(in) :: BC
    call fatalError('setBoundaryConditionsXPlane','Boundary conditions may not be set for a plane surface')
  end subroutine setBoundaryConditionsXPlane

  !!
  !! Apply boundary transformations
  !!
  subroutine boundaryTransformXPlane(self, r, u, isVacuum)
    class(xPlane), intent(in)                  :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u
    logical(defBool), intent(inout)            :: isVacuum

    if (self % isVacuum) then
      isVacuum = .TRUE.
    else if (self % isPeriodic) then
      r = r + self % periodicTranslation
    else if (self % isReflective) then
      call self % reflectiveTransform(r,u)
    else
      call fatalError('boundaryTransform, xPlane','No boundary condition applied to surface')
    end if

  end subroutine boundaryTransformXPlane

!!
!! Y-Plane procedures
!!
  subroutine initYPlane(self, y0, id, name)
    class(yPlane), intent(inout)            :: self
    real(defReal), intent(in)               :: y0
    integer(shortInt), intent(in), optional :: id
    character(*), optional, intent(in)      :: name

    self % y0 = y0
    if(present(id)) self % id = id
    if(present(name)) self % name = name

  end subroutine initYPlane

  function evaluateYPlane(self, r) result(res)
    class(yPlane), intent(in)               :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: res

    res = r(2) - self % y0

  end function evaluateYPlane

  function distanceToYPlane(self,r,u)result(distance)
    class(yPlane), intent(in)               :: self
    real(defReal), dimension(3), intent(in) :: r, &
                                               u
    real(defReal)                           :: distance, v

    v = u(2)
    distance = self%y0 - r(2)
    if ((v==ZERO) .OR. (abs(distance) < surface_tol)) then
      distance = INFINITY
      return
    end if

    distance = distance/v
    if (distance < ZERO) distance = INFINITY
    return

  end function distanceToYPlane

  subroutine reflectiveTransformYPlane(self,r,u)
    class(yPlane), intent(in)                  :: self
    real(defReal), dimension(3), intent(inout) :: r, &
                                                  u
    real(defReal) :: yDisplacement

    ! Reflect the particle y-coordinate across the plane
    yDisplacement = r(2) - self%y0
    r(2) = r(2) - TWO*yDisplacement

    ! Reflect the particle direction (independent of intersection point for plane)
    u(2) = -u(2)

  end subroutine reflectiveTransformYPlane

  function normalVectorYPlane(self,r)result(normal)
    class(yPlane), intent(in)               :: self
    real(defReal), dimension(3)             :: normal
    real(defReal), dimension(3), intent(in) :: r

    normal = [ZERO, ONE, ZERO]
    return

  end function normalVectorYPlane

  !!
  !! Give an error: this routine should not be called for a non-compound surface
  !!
  function whichYPlane(self, r, u) result(surfPointer)
    class(yPlane), intent(in)               :: self
    real(defReal), dimension(3), intent(in) :: r, u
    class(surface), pointer                 :: surfPointer
    call fatalError('whichYPlane','This function should never be called for a simple surface')
  end function whichYPlane

  !!
  !! Crash on attempting to set boundary conditions for a plane object
  !!
  subroutine setBoundaryConditionsYPlane(self, BC)
    class(yPlane), intent(inout)                :: self
    integer(shortInt), dimension(6), intent(in) :: BC
    call fatalError('setBoundaryConditionsYPlane','Boundary conditions may not be set for a plane surface')
  end subroutine setBoundaryConditionsYPlane

  !!
  !! Apply boundary transformations
  !!
  subroutine boundaryTransformYPlane(self, r, u, isVacuum)
    class(yPlane), intent(in)                  :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u
    logical(defBool), intent(inout)            :: isVacuum

    if (self % isVacuum) then
      isVacuum = .TRUE.
    else if (self % isPeriodic) then
      r = r + self % periodicTranslation
    else if (self % isReflective) then
      call self % reflectiveTransform(r,u)
    else
      call fatalError('boundaryTransform, yPlane','No boundary condition applied to surface')
    end if

  end subroutine boundaryTransformYPlane

!!
!! Z-Plane procedures
!!
  subroutine initZPlane(self, z0, id, name)
    class(zPlane), intent(inout)            :: self
    real(defReal), intent(in)               :: z0
    integer(shortInt), intent(in), optional :: id
    character(*), optional, intent(in)      :: name

    self % z0 = z0
    if(present(id)) self % id = id
    if(present(name)) self % name = name

  end subroutine initZPlane

  function evaluateZPlane(self, r) result(res)
    class(zPlane), intent(in)               :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: res

    res = r(3) - self % z0

  end function evaluateZPlane

  function distanceToZPlane(self,r,u)result(distance)
    class(zPlane), intent(in)               :: self
    real(defReal), dimension(3), intent(in) :: r, &
                                               u
    real(defReal) :: distance, w

    w = u(3)
    distance = self%z0 - r(3)
    if ((w==ZERO) .OR. (abs(distance) < surface_tol)) then
      distance = INFINITY
      return
    end if

    distance = (self % z0 - r(3))/w
    if (distance < ZERO) distance = INFINITY
    return

  end function distanceToZPlane

  subroutine reflectiveTransformZPlane(self,r,u)
    class(zPlane), intent(in)                  :: self
    real(defReal), dimension(3), intent(inout) :: r, &
                                                  u
    real(defReal) :: zDisplacement

    ! Reflect the particle z-coordinate across the plane
    zDisplacement = r(3) - self%z0
    r(3) = r(3) - TWO*zDisplacement

    ! Reflect the particle direction (independent of intersection point for plane)
    u(3) = -u(3)

  end subroutine reflectiveTransformZPlane

  function normalVectorZPlane(self,r)result(normal)
    class(zPlane), intent(in)               :: self
    real(defReal), dimension(3)             :: normal
    real(defReal), dimension(3), intent(in) :: r

    normal = [ZERO, ZERO, ONE]
    return

  end function normalVectorZPlane

  !!
  !! Give an error: this routine should not be called for a non-compound surface
  !!
  function whichZPlane(self, r, u) result(surfPointer)
    class(zPlane), intent(in)               :: self
    real(defReal), dimension(3), intent(in) :: r, u
    class(surface), pointer                 :: surfPointer
    call fatalError('whichZPlane','This function should never be called for a simple surface')
  end function whichZPlane

  !!
  !! Crash on attempting to set boundary conditions for a plane object
  !!
  subroutine setBoundaryConditionsZPlane(self, BC)
    class(zPlane), intent(inout)                :: self
    integer(shortInt), dimension(6), intent(in) :: BC
    call fatalError('setBoundaryConditionsZPlane','Boundary conditions may not be set for a plane surface')
  end subroutine setBoundaryConditionsZPlane

  !!
  !! Apply boundary transformations
  !!
  subroutine boundaryTransformZPlane(self, r, u, isVacuum)
    class(zPlane), intent(in)                  :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u
    logical(defBool), intent(inout)            :: isVacuum

    if (self % isVacuum) then
      isVacuum = .TRUE.
    else if (self % isPeriodic) then
      r = r + self % periodicTranslation
    else if (self % isReflective) then
      call self % reflectiveTransform(r,u)
    else
      call fatalError('boundaryTransform, zPlane','No boundary condition applied to surface')
    end if

  end subroutine boundaryTransformZPlane

!!
!! General Plane procedures
!!
  subroutine initPlane(self, coeff, id, name)
    class(plane), intent(inout)             :: self
    real(defReal), dimension(4), intent(in) :: coeff
    integer(shortInt), intent(in), optional :: id
    character(*), optional, intent(in)      :: name

    self % coeff = coeff
    if(present(id)) self % id = id
    if(present(name)) self % name = name

  end subroutine initPlane

  function evaluatePlane(self, r) result(res)
    class(plane), intent(in)                :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: res

    res = r(1)*self%coeff(1) + r(2)*self%coeff(2) + r(3)*self%coeff(3) - self%coeff(4)

  end function evaluatePlane

  function distanceToPlane(self,r,u) result(distance)
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

  end function distanceToPlane

  subroutine reflectiveTransformPlane(self,r,u)
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

  end subroutine reflectiveTransformPlane

  function normalVectorPlane(self,r) result(normal)
    class(plane), intent(in)                :: self
    real(defReal), dimension(3)             :: normal
    real(defReal), dimension(3), intent(in) :: r

    normal = [self%coeff(1), self%coeff(2), self%coeff(3)]
    return

  end function normalVectorPlane

  !!
  !! Give an error: this routine should not be called for a non-compound surface
  !!
  function whichPlane(self, r, u) result(surfPointer)
    class(plane), intent(in)                :: self
    real(defReal), dimension(3), intent(in) :: r, u
    class(surface), pointer                 :: surfPointer
    call fatalError('whichPlane','This function should never be called for a simple surface')
  end function whichPlane

  !!
  !! Crash on attempting to set boundary conditions for a plane object
  !!
  subroutine setBoundaryConditionsPlane(self, BC)
    class(plane), intent(inout)                 :: self
    integer(shortInt), dimension(6), intent(in) :: BC
    call fatalError('setBoundaryConditionsPlane','Boundary conditions may not be set for a plane surface')
  end subroutine setBoundaryConditionsPlane

  !!
  !! Apply boundary transformation
  !!
  subroutine boundaryTransformPlane(self, r, u, isVacuum)
    class(plane), intent(in)                   :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u
    logical(defBool), intent(inout)            :: isVacuum

    call fatalError('boundaryTransform, plane','This surface does not support boundary conditions')

  end subroutine boundaryTransformPlane

end module plane_class
