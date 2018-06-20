module zCylinder_class

  use numPrecision
  use universalVariables
  use genericProcedures, only : fatalError
  use dictionary_class,  only : dictionary
  use surface_inter,     only : surface

  implicit none
  private

  !!
  !! Constructor
  !!
  interface zCylinder
    module procedure zCylinder_fromDict
  end interface

  !!
  !! Cylinder aligned with z-axis
  !!
  type, public, extends(surface) :: zCylinder
    real(defReal)                :: rSquared = ZERO             ! squared radius of the cylinder
    real(defReal)                :: radius = ZERO               ! radius of the cylinder
    real(defReal), dimension(3)  :: origin = [ZERO, ZERO, ZERO] ! 2D origin of cylinder assumed parallel to z-axis

  contains
    procedure :: init
    procedure :: evaluate
    procedure :: distanceToSurface
    procedure :: reflectiveTransform
    procedure :: normalVector
    procedure :: whichSurface
    procedure :: setBoundaryConditions
    procedure :: boundaryTransform

  end type zCylinder


contains

  !!
  !! Initialise Z-Cylinder from components
  !!
  subroutine init(self, radius, origin, id, name)
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

  end subroutine init

  !!
  !! Returns and initialised instance of zCylinder from dictionary and name
  !!
  function zCylinder_fromDict(dict,name) result(new)
    class(dictionary), intent(in)  :: dict
    character(nameLen), intent(in) :: name
    type(zCylinder)                :: new
    integer(shortInt)              :: id
    real(defReal)                  :: radius
    real(defReal), dimension(3)    :: origin
    character(100),parameter :: Here ='zCylinder_fromDict ( zCylinder_class.f90)'

    id = dict % getInt('id')
    if(id < 1) call fatalError(Here,'Invalid surface id provided')

    radius = dict % getReal('radius')
    origin = dict % getRealArray('origin')
    call new % init(radius, origin, id, name)

  end function zCylinder_fromDict


  !!
  !! Evaluate remainder of cylinder equation
  !!
  function evaluate(self, r) result(res)
    class(zCylinder), intent(in)            :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: res

    res = (r(1) - self%origin(1))**2 + (r(2) - self%origin(2))**2 - self%rSquared

  end function evaluate

  !!
  !! Calculate distance to cylinder along direction u
  !!
  function distanceToSurface(self, r, u) result(distance)
    class(zCylinder), intent(in)            :: self
    real(defReal), dimension(3), intent(in) :: r, u
    real(defReal)                           :: xBar, yBar
    real(defReal)                           :: k, c, a
    real(defReal)                           :: discriminant, distance

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

  end function distanceToSurface

  !!
  !! Apply reflective BC
  !! ***THIS WILL NOT WORK IF USED - NEED TO ADD A DISTANCE ARGUMENT
  !!
  subroutine reflectiveTransform(self,r,u)
    class(zCylinder), intent(in)               :: self
    real(defReal), dimension(3), intent(inout) :: r, u
    character(100),parameter :: Here ='reflectiveTransform ( zCylinder_class.f90)'

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
    call fatalError(Here,'Cylinders may not have reflective boundaries')

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

  end subroutine reflectiveTransform

  !!
  !! Return normal to the cylinder
  !!
  function normalVector(self,r) result(normal)
    class(zCylinder), intent(in)            :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3)             :: normal

    normal = [TWO*(r(1) - self%origin(1)), TWO*(r(2) - self%origin(2)), ZERO]

  end function normalVector

  !!
  !! Give an error: this routine should not be called for a non-compound surface
  !!
  function whichSurface(self, r, u) result(surfPointer)
    class(zCylinder), intent(in)            :: self
    real(defReal), dimension(3), intent(in) :: r, u
    class(surface), pointer                 :: surfPointer
    character(100),parameter :: Here ='whichSurface ( zCylinder_class.f90)'

    call fatalError(Here,'This function should never be called for a simple surface')

  end function whichSurface

  !!
  !! Set boundary conditions for a zCylinder: may only be vacuum
  !!
  subroutine setBoundaryConditions(self, BC)
    class(zCylinder), intent(inout)             :: self
    integer(shortInt), dimension(6), intent(in) :: BC
    character(100),parameter :: Here ='setBoundaryConditions ( zCylinder_class.f90)'

    if (any(BC /= vacuum)) then
      call fatalError(Here,'Cylinder boundaries may only be vacuum')

    else
      self % isVacuum = .TRUE.

    end if

  end subroutine setBoundaryConditions

  !!
  !! Apply boundary transformation
  !!
  subroutine boundaryTransform(self, r, u, isVacuum)
    class(zCylinder), intent(in) :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u
    logical(defBool), intent(inout)            :: isVacuum
    character(100),parameter :: Here ='boundaryTransform ( zCylinder_class.f90)'

    if (self % isVacuum) then
      isVacuum = .TRUE.

    else
      call fatalError(Here,'This should only be called for a cylinder with vacuum boundaries')

    end if

  end subroutine boundaryTransform

    
end module zCylinder_class
