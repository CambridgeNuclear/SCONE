module yCylinder_class

  use numPrecision
  use genericProcedures, only : fatalError
  use universalVariables

  use surface_inter,     only : surface

  implicit none
  private

  !!
  !! Cylinder aligned with Y-axis
  !!
  type, public, extends(surface) :: yCylinder
    real(defReal)                :: rSquared = ZERO             ! squared radius of the cylinder
    real(defReal)                :: radius = ZERO               ! radius of the cylinder
    real(defReal), dimension(3)  :: origin = [ZERO, ZERO, ZERO] ! 2D origin of cylinder assumed parallel to y-axis

  contains
    procedure :: init
    procedure :: evaluate
    procedure :: distanceToSurface
    procedure :: reflectiveTransform
    procedure :: normalVector
    procedure :: whichSurface
    procedure :: setBoundaryConditions
    procedure :: boundaryTransform

  end type yCylinder

contains
  !!
  !! Initialise Y-Cylinder from components
  !!
  subroutine init(self, radius, origin, id, name)
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

  end subroutine init

  !!
  !! Evaluate remainder of cylinder equation
  !!
  function evaluate(self, r) result(res)
    class(yCylinder), intent(in)            :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: res

    res = (r(1) - self%origin(1))**2 + (r(3) - self%origin(3))**2 - self%rSquared
    return

  end function evaluate

  !!
  !! Calculate distance to cylinder along direction u
  !!
  function distanceToSurface(self, r, u) result(distance)
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

  end function distanceToSurface

  !!
  !! Apply reflective BC
  !! ***THIS WILL NOT WORK IF USED - NEED TO ADD A DISTANCE ARGUMENT
  !!
  subroutine reflectiveTransform(self,r,u)
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

  end subroutine reflectiveTransform

  !!
  !! Return normal to the cylinder
  !!
  function normalVector(self,r) result(normal)
    class(yCylinder), intent(in)            :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3)             :: normal

    normal = [TWO*(r(1) - self%origin(1)), ZERO, TWO*(r(3) - self%origin(3))]

  end function normalVector

  !!
  !! Give an error: this routine should not be called for a non-compound surface
  !!
  function whichSurface(self, r, u) result(surfPointer)
    class(yCylinder), intent(in)            :: self
    real(defReal), dimension(3), intent(in) :: r, u
    class(surface), pointer                 :: surfPointer

    call fatalError('whichYCylinder','This function should never be called for a simple surface')

  end function whichSurface

  !!
  !! Set boundary conditions for an yCylinder: may only be vacuum
  !!
  subroutine setBoundaryConditions(self, BC)
    class(yCylinder), intent(inout)             :: self
    integer(shortInt), dimension(6), intent(in) :: BC

    if (any(BC /= vacuum)) then
      call fatalError('setBoundaryConditionsYCylinder','Cylinder boundaries may only be vacuum')
    else
      self % isVacuum = .TRUE.
    end if
  end subroutine setBoundaryConditions

  !!
  !! Apply boundary transformation
  !!
  subroutine boundaryTransform(self, r, u, isVacuum)
    class(yCylinder), intent(in) :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u
    logical(defBool), intent(inout)            :: isVacuum

    if (self % isVacuum) then
      isVacuum = .TRUE.
    else
      call fatalError('boundaryTransform, yCylinder',&
      'This should only be called for a cylinder with vacuum boundaries')
    end if

  end subroutine boundaryTransform

    
end module yCylinder_class
