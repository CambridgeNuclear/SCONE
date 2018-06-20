module xCylinder_class

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
  interface xCylinder
    module procedure xCylinder_fromDict
  end interface

  !!
  !! Cylinder aligned with x-axis
  !!
  type, public, extends(surface) :: xCylinder
    real(defReal)                :: rSquared = ZERO             ! squared radius of the cylinder
    real(defReal)                :: radius = ZERO               ! radius of the cylinder
    real(defReal), dimension(3)  :: origin = [ZERO, ZERO, ZERO] ! 2D origin of cylinder assumed parallel to x-axis

  contains
    procedure :: init
    procedure :: evaluate
    procedure :: distanceToSurface
    procedure :: reflectiveTransform
    procedure :: normalVector
    procedure :: whichSurface
    procedure :: setBoundaryConditions
    procedure :: boundaryTransform

  end type xCylinder

contains

  !!
  !! Initialise X-Cylinder from components
  !!
  subroutine init(self, radius, origin, id, name)
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

  end subroutine init

  !!
  !! Returns and initialised instance of xCylinder from dictionary and name
  !!
  function xCylinder_fromDict(dict,name) result(new)
    class(dictionary), intent(in)  :: dict
    character(nameLen), intent(in) :: name
    type(xCylinder)                :: new
    integer(shortInt)              :: id
    real(defReal)                  :: radius
    real(defReal), dimension(3)    :: origin
    character(100),parameter :: Here ='xCylinder_fromDict ( xCylinder_class.f90)'

    id = dict % getInt('id')
    if(id < 1) call fatalError(Here,'Invalid surface id provided')

    radius = dict % getReal('radius')
    origin = dict % getRealArray('origin')
    call new % init(radius, origin, id, name)

  end function xCylinder_fromDict

  !!
  !! Evaluate remainder of cylinder equation
  !!
  function evaluate(self, r) result(res)
    class(xCylinder), intent(in) :: self
    real(defReal), dimension(3), intent(in):: r
    real(defReal) :: res

    res = (r(2) - self%origin(2))**2 + (r(3) - self%origin(3))**2 - self%rSquared

  end function evaluate

  !!
  !! Calculate distance to cylinder along direction u
  !!
  function distanceToSurface(self, r, u) result(distance)
    class(xCylinder), intent(in)            :: self
    real(defReal), dimension(3), intent(in) :: r, u
    real(defReal)                           :: yBar, zBar
    real(defReal)                           :: k, c, a
    real(defReal)                           :: discriminant,  distance

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

  end function distanceToSurface

  !!
  !! Apply reflective BC
  !! ***THIS WILL NOT WORK IF USED - NEED TO ADD A DISTANCE ARGUMENT
  !!
  subroutine reflectiveTransform(self,r,u)
    class(xCylinder), intent(in)               :: self
    real(defReal), dimension(3), intent(inout) :: r, u
    character(100),parameter :: Here ='reflectiveTransform ( xCylinder_class.f90)'
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
    call fatalError(Here,'Cylinders may not have reflective boundaries')

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

  end subroutine reflectiveTransform

  !!
  !! Return normal to the cylinder
  !!
  function normalVector(self,r) result(normal)
    class(xCylinder), intent(in)            :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3)             :: normal

    normal = [ZERO, TWO*(r(2) - self%origin(2)), TWO*(r(3) - self%origin(3))]

  end function normalVector

  !!
  !! Give an error: this routine should not be called for a non-compound surface
  !!
  function whichSurface(self, r, u) result(surfPointer)
    class(xCylinder), intent(in)            :: self
    real(defReal), dimension(3), intent(in) :: r, u
    class(surface), pointer                 :: surfPointer
    character(100),parameter   :: Here ='whichSurface ( xCylinder_class.f90)'

    call fatalError(Here,'This function should never be called for a simple surface')

  end function whichSurface

  !!
  !! Set boundary conditions for an xCylinder: may only be vacuum
  !!
  subroutine setBoundaryConditions(self, BC)
    class(xCylinder), intent(inout)             :: self
    integer(shortInt), dimension(6), intent(in) :: BC
    character(100),parameter :: Here ='setBoundaryConditions ( xCylinder_class.f90)'

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
    class(xCylinder), intent(in) :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u
    logical(defBool), intent(inout)            :: isVacuum
    character(100),parameter :: Here ='boundaryTransform ( xCylinder_class.f90)'

    if (self % isVacuum) then
      isVacuum = .TRUE.

    else
      call fatalError(Here,'This should only be called for a cylinder with vacuum boundaries')

    end if

  end subroutine boundaryTransform

    
end module xCylinder_class
