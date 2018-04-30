module sphere_class      

  use numPrecision
  use genericProcedures
  use universalVariables

  use surface_class
  use plane_class

  implicit none
  private

  !
  ! Sphere
  ! f = (x-x0)^2 + (y-y0)^2 + (z-z0)^2 = R^2
  !
  type, public, extends (surface) :: sphere
    private
    real(defReal), dimension(3) :: origin
    real(defReal) :: rSquared = 0.0
    real(defReal) :: radius = 0.0
  contains
    procedure :: init => initSphere         ! provide sphere with origin and radius
    procedure :: evaluate => evaluateSphere
    procedure :: distanceToSurface => distanceToSphere
    procedure :: reflectiveTransform => reflectiveTransformSphere
    !procedure :: periodicBoundary
    procedure :: normalVector => normalVectorSphere
    procedure :: whichSurface
    procedure :: setBoundaryConditions => setBoundaryConditionsSphere
  end type sphere

contains

  !
  ! Given an origin and radius, create a sphere object
  !
  subroutine initSphere(self, origin, radius, id, name)
    implicit none
    class(sphere), intent(inout) :: self
    real(defReal), dimension(3), intent(in) :: origin
    real(defReal), intent (in) :: radius
    integer(shortInt), intent(in), optional :: id
    character(*), optional, intent(in) :: name

    self % origin = origin
    self % radius = radius
    self % rSquared = radius*radius
    if(present(id)) self % id = id
    if(present(name)) self % name = name

  end subroutine initSphere

  !
  ! Check which sphere halfspace a point occupies - returns -1, 0, or 1
  !
  ! Calculate squared difference between point and origin of sphere
  ! Insert values into sphere equation
  !
  function evaluateSphere(self, r) result(res)
    implicit none
    class(sphere), intent(in) :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3) :: diff
    real(defReal) :: res

    diff=(r - self % origin)**2
    res = diff(1) + diff(2) + diff(3) - self % rSquared

  end function evaluateSphere

  !
  ! Calculate the neutron's distance to the sphere surface
  ! Solves quadratic equation: d^2 + 2kd + c = 0 where:
  ! k = (x-x0)u + (y-y0)v + (z-z0)w
  ! c = (x-x0)^2 + (y-y0)^2 + (z-z0)^2 - R^2
  !
  function distanceToSphere(self,r,u)result(distance)
    implicit none
    class(sphere), intent(in) :: self
    real(defReal), dimension(3), intent(in) :: r, &
                                               u
    real(defReal), dimension(3) :: rBar
    real(defReal) :: k, &
                     c, &
                     discriminant, &
                     distance

    rBar = r - self%origin
    k = dotProduct(rBar,u)
    c = dotProduct(rBar,rBar) - self%rSquared
    discriminant = k*k - c

    ! Particle does not intersect the surface
    if(discriminant < 0.0 ) then
      distance = INFINITY
      return
    ! Particle is on the sphere - the distance will be either positive or negative
    else if (abs(c) < surface_tol) then
      if (k >= 0.0) then
        distance = INFINITY
        return
      else
        distance = -k + sqrt(discriminant)
        return
      end if
    ! Particle is inside the surface
    else if (c < 0.0) then
      distance = -k + sqrt(discriminant)
      return
    ! The neutron is outside the surface
    else
      distance = -k - sqrt(discriminant)
      if (distance < 0.0) distance = INFINITY
      return
    end if

  end function distanceToSphere

  !
  ! Perform a co-ordinate transform on a particle to apply reflective boundary condition
  ! Used as the standard reflection is not particularly agreeable when used with delta tracking
  !
  ! Calculate particle's initial angle to the centre point of the sphere
  ! Then find the particle's new position after being moved by a given distance
  ! Calculate the outward normal where the particle crossed the surface
  ! Reflect the particle in the plane given by this outward normal
  !
  ! THIS WILL NOT WORK IF USED - NEED TO ADD A DISTANCE ARGUMENT
  subroutine reflectiveTransformSphere(self, r, u)
    implicit none
    class(sphere), intent(in) :: self
    real(defReal), dimension(3), intent(inout) :: r, &
                                                  u
    real(defReal), dimension(3) :: Ovector, &    ! vector towards origin from neutron
                                   normal, &     ! normal vector of the reflection plane
                                   perpVector, & ! normal vector to the plane in which motion/rotation is occurring
                                   intersect     ! point at which the particle intersects the sphere
    real(defReal) :: cosGamma, &         ! cosine of the half angle of reflection
                     perpDistanceOrig, & ! perpendicular distance between neutron path and origin
                     sinGamma, &         ! sine of the half angle of reflection
                     perpDistance, &     ! distance of a point outside surface to reflective plane
                     dOrigin             ! distance from the starting point to the origin

    ! Reflective transforms will not be allowed to occur in geometries other than planes
    call fatalError('reflectiveTransform','Spheres may not have reflective boundaries')

    ! Calculate the normalised origin vector
    Ovector= self % origin - r
    dOrigin=norm2(Ovector)
    Ovector = Ovector / dOrigin

    ! Obtain a perpendicular vector by taking the cross product
    ! Rotation will take place around this axis
    perpVector = crossProduct(Ovector,u)

    ! Calculate sinGamma using the cross-product of the two vectors and the sine law
    sinGamma = (dOrigin/self%radius)*norm2(perpVector)
    cosGamma = sqrt(1 - sinGamma*sinGamma)

    ! Rotate the particle direction by gamma to obtain the reflection plane normal
    ! Use Rodrigues' rotation formula
    ! (Neglect [1 - cosGamma] term in formula due to taking a dotProduct of perp. vectors)
    normal = cosGamma*u + sinGamma*crossProduct(perpVector,u)

    ! Obtain the point at which the intersection occurs from the normal, radius and origin
    intersect = self%origin + normal * self%radius

    ! Calculate particle position outside the sphere
    !r = r + u*distance
    !MAY NEED TO FIX

    ! Calculate the perpendicular distance to the reflecting plane
    perpDistance = abs(dotProduct(normal,r-intersect))

    ! Translate the particle position across the plane
    r = r - 2*perpDistance*normal

    ! Reflect the particle direction (independent of intersection point for plane)
    u = u - 2*dotProduct(normal,u)*normal

  end subroutine reflectiveTransformSphere

  !
  ! Supply the normal vector given a point on the sphere
  !
  function normalVectorSphere(self, r) result(normal)
    implicit none
    class(sphere), intent(in) :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3) :: normal

    normal = 2.0*(r - self % origin)

  end function normalVectorSphere

  !
  ! Give an error: this routine should not be called for a non-compound surface
  !
  function whichSurface(self, r, u) result(surfPointer)
    implicit none
    class(sphere), intent(in) :: self
    real(defReal), dimension(3), intent(in) :: r, u
    class(surface), pointer :: surfPointer
    call fatalError('whichSurface (sphere)','This function should never be called for a simple surface')
  end function whichSurface

  !!
  !! Set boundary conditions for a sphere: may only be vacuum
  !!
  subroutine setBoundaryConditionsSphere(self, BC)
    class(sphere), intent(inout) :: self
    integer(shortInt), dimension(6), intent(in) :: BC

    if (any(BC /= vacuum)) then
      call fatalError('setBoundaryConditionsSphere','Sphere boundaries may only be vacuum')
    else
      self % isVacuum = .TRUE.
    end if
  end subroutine setBoundaryConditionsSphere

end module sphere_class
