!
! Module containing cylinders parallel to the x, y, and z axis
!
module cylinder_class
  use numPrecision
  use genericProcedures, only : fatalError
  use universalVariables

  use surface_class

  implicit none
  private

  ! X-cylinder definition
  type, public, extends(surface) :: xCylinder
    real(defReal)                :: rSquared = ZERO        ! squared radius of the cylinder
    real(defReal)                :: radius = ZERO          ! radius of the cylinder
    real(defReal), dimension(3)  :: origin = [ZERO, ZERO, ZERO] ! 2D origin of cylinder assumed parallel to x-axis
  contains
    procedure :: init => initXCyl
    procedure :: evaluate => evaluateXCyl
    procedure :: distanceToSurface => distanceToXCyl
    procedure :: reflectiveTransform => reflectiveTransformXCyl
    procedure :: normalVector => normalVectorXCyl
    procedure :: whichSurface => whichXCylinder
    procedure :: setBoundaryConditions => setBoundaryConditionsXCylinder
  end type xCylinder

  ! Y-cylinder definition
  type, public, extends(surface) :: yCylinder
    real(defReal)                :: rSquared = ZERO             ! squared radius of the cylinder
    real(defReal)                :: radius = ZERO               ! radius of the cylinder
    real(defReal), dimension(3)  :: origin = [ZERO, ZERO, ZERO] ! 2D origin of cylinder assumed parallel to y-axis
  contains
    procedure :: init => initYCyl
    procedure :: evaluate => evaluateYCyl
    procedure :: distanceToSurface => distanceToYCyl
    procedure :: reflectiveTransform => reflectiveTransformYCyl
    procedure :: normalVector => normalVectorYCyl
    procedure :: whichSurface => whichYCylinder
    procedure :: setBoundaryConditions => setBoundaryConditionsYCylinder
  end type yCylinder

  ! Z-cylinder definition
  type, public, extends(surface) :: zCylinder
    real(defReal)                :: rSquared = ZERO             ! squared radius of the cylinder
    real(defReal)                :: radius = ZERO               ! radius of the cylinder
    real(defReal), dimension(3)  :: origin = [ZERO, ZERO, ZERO] ! 2D origin of cylinder assumed parallel to z-axis
  contains
    procedure :: init => initZCyl
    procedure :: evaluate => evaluateZCyl
    procedure :: distanceToSurface => distanceToZCyl
    procedure :: reflectiveTransform => reflectiveTransformZCyl
    procedure :: normalVector => normalVectorZCyl
    procedure :: whichSurface => whichZCylinder
    procedure :: setBoundaryConditions => setBoundaryConditionsZCylinder
  end type zCylinder

contains

!
! X-cylinder procedures
!

  subroutine initXCyl(self, radius, origin, id, name)
    class(xCylinder), intent(inout)         :: self
    real(defReal), intent(in)               :: radius
    real(defReal), dimension(3), intent(in) :: origin
    integer(shortInt), intent(in), optional :: id
    character(*), optional, intent(in)      :: name

    self % radius = radius
    if(radius < surface_tol) &
    call fatalError('init, xCylinder','Radius must be greater than surface tolerance')
    self % rSquared = radius*radius
    self % origin = origin
    if(present(id)) self % id = id
    if(present(name)) self % name = name

  end subroutine initXCyl

  function evaluateXCyl(self, r) result(res)
    class(xCylinder), intent(in) :: self
    real(defReal), dimension(3), intent(in):: r
    real(defReal) :: res

    res = (r(2) - self%origin(2))**2 + (r(3) - self%origin(3))**2 - self%rSquared
    return

  end function evaluateXCyl

  function distanceToXCyl(self, r, u) result(distance)
    class(xCylinder), intent(in)            :: self
    real(defReal), dimension(3), intent(in) :: r, &
                                               u
    real(defReal)                           :: yBar, &
                                               zBar, &
                                               k, &
                                               c, &
                                               a, &
                                               discriminant, &
                                               distance

    yBar = r(2) - self%origin(2)
    zBar = r(3) - self%origin(3)

    k = yBar*u(2) + zBar*u(3)
    a = ONE - u(1)*u(1) ! = u(2)*u(2) + u(3)*u(3)
    c = yBar*yBar + zBar*zBar - self%rSquared
    discriminant = k*k - a*c

    ! No intersection
    if ((a==ZERO).OR.(discriminant<ZERO)) then
      distance = INFINITY
      return
    ! Particle is on the cylinder - the distance will be either positive or negative
    else if (abs(c) < surface_tol) then
      if (k >= ZERO) then
        distance = INFINITY
        return
      else
        distance = (-k + sqrt(discriminant))/a
        return
      end if
    ! Particle is inside - take solution with + before sqrt
    else if (c < ZERO) then
      distance = (-k + sqrt(discriminant))/a
      return
    ! Particle is outside - both distances are either positive or negative
    else
      distance = (-k -sqrt(discriminant))/a
      if (distance < ZERO) distance = INFINITY
      return
    end if

  end function distanceToXCyl

  ! THIS WILL NOT WORK IF USED - NEED TO ADD A DISTANCE ARGUMENT
  subroutine reflectiveTransformXCyl(self,r,u)
    class(xCylinder), intent(in)               :: self
    real(defReal), dimension(3), intent(inout) :: r, &
                                                  u
    !real(defReal), dimension(3) :: normal, &
    !                               Ovector, &
    !                               yzVector, &
    !                               intersect
    !real(defReal) :: magSquared, &
    !                 radius, &
    !                 dOrigin, &
    !                 perpDistance, &
    !                 cosGamma, &
    !                 sinGamma

    ! Reflective transforms will not be allowed to occur in geometries other than planes
    call fatalError('reflectiveTransformXCyl, xCylinder','Cylinders may not have reflective boundaries')

    ! Construct unit vector from origin to starting point
    ! Zero the x-entry as there is no true x-origin
    !Ovector = self%origin - r
    !Ovector(1) = 0.0
    ! Calculate the distance to the origin
    !dOrigin = norm2(Ovector)
    !Ovector = Ovector/dOrigin

    ! Determine the direction vector of the particle in the yz-plane
    !yzVector = u
    !yzVector(1) = 0.0
    !yzVector = yzVector/norm2(yzVector)

    ! Calculate sinGamma using the cross-product of the two vectors and the sine law
    !sinGamma = (dOrigin/self%radius)*norm2(crossProduct(Ovector,yzVector))
    !cosGamma = sqrt(1 - sinGamma*sinGamma)

    ! Rotate the particle direction by gamma to obtain the normal vector
    !normal = yzVector
    !normal(2) = yzVector(2) * cosGamma - yzVector(3) * sinGamma
    !normal(3) = yzVector(2) * sinGamma + yzVector(3) * cosGamma

    ! Obtain the point at which the intersection occurs from the normal, radius and origin
    !intersect = self%origin + normal * self%radius

    ! Calculate particle position outside the cylinder
    !r = r + u*distance

    ! Calculate the perpendicular distance to the plane
    !perpDistance = abs(dotProduct(normal,r-intersect))

    ! Translate the particle position across the plane
    !r = r - 2*perpDistance*normal

    ! Reflect the particle direction (independent of intersection point for plane)
    !u = u - 2*dotProduct(normal,u)*normal

  end subroutine reflectiveTransformXCyl

  function normalVectorXCyl(self,r) result(normal)
    class(xCylinder), intent(in)            :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3)             :: normal

    normal = [ZERO, TWO*(r(2) - self%origin(2)), TWO*(r(3) - self%origin(3))]
    return

  end function normalVectorXCyl

  !
  ! Give an error: this routine should not be called for a non-compound surface
  !
  function whichXCylinder(self, r, u) result(surfPointer)
    class(xCylinder), intent(in)            :: self
    real(defReal), dimension(3), intent(in) :: r, u
    class(surface), pointer                 :: surfPointer
    call fatalError('whichXCylinder, xCylinder','This function should never be called for a simple surface')
  end function whichXCylinder

  !!
  !! Set boundary conditions for an xCylinder: may only be vacuum
  !!
  subroutine setBoundaryConditionsXCylinder(self, BC)
    class(xCylinder), intent(inout)             :: self
    integer(shortInt), dimension(6), intent(in) :: BC

    if (any(BC /= vacuum)) then
      call fatalError('setBoundaryConditionsXCylinder','Cylinder boundaries may only be vacuum')
    else
      self % isVacuum = .TRUE.
    end if
  end subroutine setBoundaryConditionsXCylinder

!
! Y-cylinder procedures
!

  subroutine initYCyl(self, radius, origin, id, name)
    class(yCylinder), intent(inout)         :: self
    real(defReal), intent(in)               :: radius
    real(defReal), dimension(3), intent(in) :: origin
    integer(shortInt), intent(in), optional :: id
    character(*), optional, intent(in)      :: name

    self % radius = radius
    if (radius < surface_tol) &
    call fatalError('init, yCylinder','Radius must be greater than surface tolerance')
    self % rSquared = radius*radius
    self % origin = origin
    if(present(id)) self % id = id
    if(present(name)) self % name = name

  end subroutine initYCyl

  function evaluateYCyl(self, r) result(res)
    class(yCylinder), intent(in)            :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: res

    res = (r(1) - self%origin(1))**2 + (r(3) - self%origin(3))**2 - self%rSquared
    return

  end function evaluateYCyl

  function distanceToYCyl(self, r, u) result(distance)
    class(yCylinder), intent(in)            :: self
    real(defReal), dimension(3), intent(in) :: r, &
                                               u
    real(defReal)                           :: xBar, &
                                               zBar, &
                                               k, &
                                               c, &
                                               a, &
                                               discriminant, &
                                               distance

    xBar = r(1) - self%origin(1)
    zBar = r(3) - self%origin(3)

    k = xBar*u(1) + zBar*u(3)
    a = ONE - u(2)*u(2) ! = u(1)*u(1) + u(3)*u(3)
    c = xBar*xBar + zBar*zBar - self%rSquared
    discriminant = k*k - a*c

    ! No intersection
    if ((a==ZERO).OR.(discriminant<ZERO)) then
      distance = INFINITY
      return
    ! Particle is on the cylinder - the distance will be either positive or negative
    else if (abs(c) < surface_tol) then
      if (k >= ZERO) then
        distance = INFINITY
        return
      else
        distance = (-k + sqrt(discriminant))/a
        return
      end if
    ! Particle is inside - take solution with + before sqrt
    else if (c<ZERO) then
      distance = (-k + sqrt(discriminant))/a
      return
    ! Particle is outside - both distances are either positive or negative
    else
      distance = (-k -sqrt(discriminant))/a
      if (distance<ZERO) distance = INFINITY
      return
    end if

  end function distanceToYCyl

  ! THIS WILL NOT WORK IF USED - NEED TO ADD A DISTANCE ARGUMENT
  subroutine reflectiveTransformYCyl(self,r,u)
    class(yCylinder), intent(in)               :: self
    real(defReal), dimension(3), intent(inout) :: r, &
                                                  u
    !real(defReal), dimension(3) :: normal, &
    !                               Ovector, &
    !                               xzVector, &
    !                               intersect
    !real(defReal) :: magSquared, &
    !                 radius, &
    !                 dOrigin, &
    !                 perpDistance, &
    !                 cosGamma, &
    !                 sinGamma

    ! Reflective transforms will not be allowed to occur in geometries other than planes
    call fatalError('reflectiveTransformYCyl, yCylinder','Cylinders may not have reflective boundaries')

    ! Construct unit vector from origin to starting point
    ! Zero the y-entry as there is no true y-origin
    !Ovector = self%origin - r
    !Ovector(2) = 0.0
    ! Calculate the distance to the origin
    !dOrigin = norm2(Ovector)
    !Ovector = Ovector/dOrigin

    ! Determine the direction vector of the particle in the xz-plane
    !xzVector = u
    !xzVector(2) = 0.0
    !xzVector = xzVector/norm2(xzVector)

    ! Calculate sinGamma using the cross-product of the two vectors and the sine law
    !sinGamma = (dOrigin/self%radius)*norm2(crossProduct(Ovector,xzVector))
    !cosGamma = sqrt(1 - sinGamma*sinGamma)

    ! Rotate the particle direction by gamma to obtain the normal vector
    !normal = xzVector
    !normal(1) = xzVector(1) * cosGamma - xzVector(3) * sinGamma
    !normal(3) = xzVector(1) * sinGamma + xzVector(3) * cosGamma

    ! Obtain the point at which the intersection occurs from the normal, radius and origin
    !intersect = self%origin + normal * self%radius

    ! Calculate particle position outside the cylinder
    !r = r + u*distance

    ! Calculate the perpendicular distance to the plane
    !perpDistance = abs(dotProduct(normal,r-intersect))

    ! Translate the particle position across the plane
    !r = r - 2*perpDistance*normal

    ! Reflect the particle direction (independent of intersection point for plane)
    !u = u - 2*dotProduct(normal,u)*normal

  end subroutine reflectiveTransformYCyl

  function normalVectorYCyl(self,r) result(normal)
    class(yCylinder), intent(in)            :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3)             :: normal

    normal = [TWO*(r(1) - self%origin(1)), ZERO, TWO*(r(3) - self%origin(3))]
    return

  end function normalVectorYCyl

  !
  ! Give an error: this routine should not be called for a non-compound surface
  !
  function whichYCylinder(self, r, u) result(surfPointer)
    class(yCylinder), intent(in)            :: self
    real(defReal), dimension(3), intent(in) :: r, u
    class(surface), pointer                 :: surfPointer
    call fatalError('whichYCylinder','This function should never be called for a simple surface')
  end function whichYCylinder

  !!
  !! Set boundary conditions for an yCylinder: may only be vacuum
  !!
  subroutine setBoundaryConditionsYCylinder(self, BC)
    class(yCylinder), intent(inout)             :: self
    integer(shortInt), dimension(6), intent(in) :: BC

    if (any(BC /= vacuum)) then
      call fatalError('setBoundaryConditionsYCylinder','Cylinder boundaries may only be vacuum')
    else
      self % isVacuum = .TRUE.
    end if
  end subroutine setBoundaryConditionsYCylinder

!
! Z-cylinder procedures
!

  subroutine initZCyl(self, radius, origin, id, name)
    class(zCylinder), intent(inout)         :: self
    real(defReal), intent(in)               :: radius
    real(defReal), dimension(3), intent(in) :: origin
    integer(shortInt), intent(in), optional :: id
    character(*), optional, intent(in)      :: name

    self % radius = radius
    if (radius < surface_tol) &
    call fatalError('init, zCylinder','Radius must be greater than surface tolerance')
    self % rSquared = radius*radius
    self % origin = origin
    if(present(id)) self % id = id
    if(present(name)) self % name = name

  end subroutine initZCyl

  function evaluateZCyl(self, r) result(res)
    class(zCylinder), intent(in)            :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: res

    res = (r(1) - self%origin(1))**2 + (r(2) - self%origin(2))**2 - self%rSquared
    return

  end function evaluateZCyl

  function distanceToZCyl(self, r, u) result(distance)
    class(zCylinder), intent(in)            :: self
    real(defReal), dimension(3), intent(in) :: r, &
                                               u
    real(defReal)                           :: xBar, &
                                               yBar, &
                                               k, &
                                               c, &
                                               a, &
                                               discriminant, &
                                               distance

    xBar = r(1) - self%origin(1)
    yBar = r(2) - self%origin(2)

    k = xBar*u(1) + yBar*u(2)
    a = ONE - u(3)*u(3) ! = u(1)*u(1) + u(2)*u(2)
    c = xBar*xBar + yBar*yBar - self%rSquared
    discriminant = k*k - a*c

    ! No intersection
    if ((a==ZERO).OR.(discriminant<ZERO)) then
      distance = INFINITY
      return
    ! Particle is on the cylinder - the distance will be either positive or negative
    else if (abs(c) < surface_tol) then
      if (k >= ZERO) then
        distance = INFINITY
        return
      else
        distance = (-k + sqrt(discriminant))/a
        return
      end if
    ! Particle is inside - take solution with + before sqrt
    else if (c<ZERO) then
      distance = (-k + sqrt(discriminant))/a
      return
    ! Particle is outside - both distances are either positive or negative
    else
      distance = (-k -sqrt(discriminant))/a
      if (distance<ZERO) distance = INFINITY
      return
    end if

  end function distanceToZCyl

  ! THIS WILL NOT WORK IF USED - NEED TO ADD A DISTANCE ARGUMENT
  subroutine reflectiveTransformZCyl(self,r,u)
    class(zCylinder), intent(in)               :: self
    real(defReal), dimension(3), intent(inout) :: r, &
                                                  u
    !real(defReal), dimension(3) :: normal, &
    !                               Ovector, &
    !                               xyVector, &
    !                               intersect
    !real(defReal) :: magSquared, &
    !                 radius, &
    !                 dOrigin, &
    !                 perpDistance, &
    !                 cosGamma, &
    !                 sinGamma

    ! Reflective transforms will not be allowed to occur in geometries other than planes
    call fatalError('reflectiveTransformZCyl,zCylinder','Cylinders may not have reflective boundaries')

    ! Construct unit vector from origin to starting point
    ! Zero the z-entry as there is no true z-origin
    !Ovector = self%origin - r
    !Ovector(3) = 0.0
    ! Calculate the distance to the origin
    !dOrigin = norm2(Ovector)
    !Ovector = Ovector/dOrigin

    ! Determine the direction vector of the particle in the xy-plane
    !xyVector = u
    !xyVector(3) = 0.0
    !xyVector = xyVector/norm2(xyVector)

    ! Calculate sinGamma using the cross-product of the two vectors and the sine law
    !sinGamma = (dOrigin/self%radius)*norm2(crossProduct(Ovector,xyVector))
    !cosGamma = sqrt(1 - sinGamma*sinGamma)

    ! Rotate the particle direction by gamma to obtain the normal vector
    !normal = xyVector
    !normal(1) = xyVector(1) * cosGamma - xyVector(2) * sinGamma
    !normal(2) = xyVector(1) * sinGamma + xyVector(2) * cosGamma

    ! Obtain the point at which the intersection occurs from the normal, radius and origin
    !intersect = self%origin + normal * self%radius

    ! Calculate particle position outside the cylinder
    !r = r + u*distance

    ! Calculate the perpendicular distance to the plane
    !perpDistance = abs(dotProduct(normal,r-intersect))

    ! Translate the particle position across the plane
    !r = r - 2*perpDistance*normal

    ! Reflect the particle direction (independent of intersection point for plane)
    !u = u - 2*dotProduct(normal,u)*normal

  end subroutine reflectiveTransformZCyl

  function normalVectorZCyl(self,r) result(normal)
    class(zCylinder), intent(in)            :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3)             :: normal

    normal = [TWO*(r(1) - self%origin(1)), TWO*(r(2) - self%origin(2)), ZERO]
    return

  end function normalVectorZCyl

  !
  ! Give an error: this routine should not be called for a non-compound surface
  !
  function whichZCylinder(self, r, u) result(surfPointer)
    class(zCylinder), intent(in)            :: self
    real(defReal), dimension(3), intent(in) :: r, u
    class(surface), pointer                 :: surfPointer
    call fatalError('whichZCylinder, zCylinder','This function should never be called for a simple surface')
  end function whichZCylinder

  !!
  !! Set boundary conditions for a zCylinder: may only be vacuum
  !!
  subroutine setBoundaryConditionsZCylinder(self, BC)
    class(zCylinder), intent(inout)             :: self
    integer(shortInt), dimension(6), intent(in) :: BC

    if (any(BC /= vacuum)) then
      call fatalError('setBoundaryConditionsZCylinder','Cylinder boundaries may only be vacuum')
    else
      self % isVacuum = .TRUE.
    end if
  end subroutine setBoundaryConditionsZCylinder

end module cylinder_class
