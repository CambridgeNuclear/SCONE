module sphere_class      

  use numPrecision
  use universalVariables
  use genericProcedures, only : fatalError, dotProduct
  use dictionary_class,  only : dictionary
  use surface_inter,     only : surface, printSurfDef

  implicit none
  private

  !!
  !! Constants describing surface properties
  !!
  character(nameLen),parameter :: TYPE_NAME    = 'sphere'
  logical(defBool),parameter   :: ACCEPTS_BC   = .true.

  !!
  !! Constructor
  !!
  interface sphere
    module procedure sphere_fromDict
  end interface

  !!
  !! Sphere
  !! f = (x-x0)^2 + (y-y0)^2 + (z-z0)^2 = R^2
  !!
  type, public, extends (surface) :: sphere
    private
    real(defReal), dimension(3) :: origin
    real(defReal) :: rSquared = ZERO
    real(defReal) :: radius = ZERO

  contains
    procedure :: init
    procedure :: evaluate
    procedure :: type
    procedure :: getDef
    procedure :: cannotBeBoundary
    procedure :: distanceToSurface
    procedure :: reflectiveTransform
    procedure :: normalVector
    procedure :: whichSurface
    procedure :: setBoundaryConditions
    procedure :: boundaryTransform

  end type sphere

contains

  !!
  !! Given an origin and radius, create a sphere object
  !!
  subroutine init(self, origin, radius, id, name)
    class(sphere), intent(inout)            :: self
    real(defReal), dimension(3), intent(in) :: origin
    real(defReal), intent (in)              :: radius
    integer(shortInt), intent(in), optional :: id
    character(*), optional, intent(in)      :: name

    self % origin = origin
    self % radius = radius
    if (radius < surface_tol) &
    call fatalError('initSphere, sphere','Radius must be greater than surface tolerance')
    self % rSquared = radius*radius
    if(present(id)) self % id = id
    if(present(name)) self % name = name

    ! Hash and store surface definition
    call self % hashSurfDef()

  end subroutine init

  !!
  !! Returns and initialised instance of sphere from dictionary and name
  !!
  function sphere_fromDict(dict,name) result(new)
    class(dictionary), intent(in)  :: dict
    character(nameLen), intent(in) :: name
    type(sphere)                   :: new
    integer(shortInt)              :: id
    real(defReal)                  :: radius
    real(defReal), dimension(3)    :: origin
    character(100), parameter :: Here = 'sphere_fromDict (sphere_class.f90)'

    id = dict % getInt('id')
    if(id < 1) call fatalError(Here,'Invalid surface id provided')

    radius = dict % getReal('radius')
    origin = dict % getRealArray('origin')

    call new % init(origin, radius, id, name)

  end function sphere_fromDict

  !!
  !! Calculate squared difference between point and origin of sphere
  !! Insert values into sphere equation
  !!
  function evaluate(self, r) result(res)
    class(sphere), intent(in)               :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3)             :: diff
    real(defReal)                           :: res

    diff=(r - self % origin)**2
    res = diff(1) + diff(2) + diff(3) - self % rSquared

  end function evaluate

  !!
  !! Return parameter character containing TYPE NAME
  !!
  function type(self)
    class(sphere), intent(in) :: self
    character(nameLen)        :: type

    type = TYPE_NAME

  end function type

  !!
  !! Return string with definition of this surface
  !!
  pure subroutine getDef(self,string)
    class(sphere), intent(in)              :: self
    character(:),allocatable,intent(inout) :: string

    string = printSurfDef(TYPE_NAME, [self % radius, self % origin])

  end subroutine getDef

  !!
  !! Override base type function to returns .false.
  !!
  function cannotBeBoundary(self) result(itCant)
    class(sphere), intent(in) :: self
    logical(defBool)          :: itCant

    itCant = .not.ACCEPTS_BC

  end function cannotBeBoundary

  !!
  !! Calculate the neutron's distance to the sphere surface
  !! Solves quadratic equation: d^2 + 2kd + c = 0 where:
  !! k = (x-x0)u + (y-y0)v + (z-z0)w
  !! c = (x-x0)^2 + (y-y0)^2 + (z-z0)^2 - R^2
  !!
  function distanceToSurface(self,r,u)result(distance)
    class(sphere), intent(in)               :: self
    real(defReal), dimension(3), intent(in) :: r, &
                                               u
    real(defReal), dimension(3)             :: rBar
    real(defReal)                           :: k, &
                                               c, &
                                               discriminant, &
                                               distance

    rBar = r - self%origin
    k = dotProduct(rBar,u)
    c = dotProduct(rBar,rBar) - self%rSquared
    discriminant = k*k - c

    ! Particle does not intersect the surface
    if(discriminant < ZERO ) then
      distance = INFINITY
      return
    ! Particle is on the sphere - the distance will be either positive or negative
    else if (abs(c) < surface_tol) then
      if (k >= ZERO) then
        distance = INFINITY
        return
      else
        distance = -k + sqrt(discriminant)
        return
      end if
    ! Particle is inside the surface
    else if (c < ZERO) then
      distance = -k + sqrt(discriminant)
      return
    ! The neutron is outside the surface
    else
      distance = -k - sqrt(discriminant)
      if (distance < ZERO) distance = INFINITY
      return
    end if

  end function distanceToSurface

  !!
  !! Perform a co-ordinate transform on a particle to apply reflective boundary condition
  !! Used as the standard reflection is not particularly agreeable when used with delta tracking
  !!
  !! Calculate particle's initial angle to the centre point of the sphere
  !! Then find the particle's new position after being moved by a given distance
  !! Calculate the outward normal where the particle crossed the surface
  !! Reflect the particle in the plane given by this outward normal
  !!
  !! THIS WILL NOT WORK IF USED - NEED TO ADD A DISTANCE ARGUMENT
  subroutine reflectiveTransform(self, r, u)
    class(sphere), intent(in)                  :: self
    real(defReal), dimension(3), intent(inout) :: r, u
    character(100), parameter :: Here = 'reflectiveTransform (sphere_class.f90)'
    !real(defReal), dimension(3) :: Ovector, &    ! vector towards origin from neutron
    !                               normal, &     ! normal vector of the reflection plane
    !                               perpVector, & ! normal vector to the plane in which motion/rotation is occurring
    !                               intersect     ! point at which the particle intersects the sphere
    !real(defReal) :: cosGamma, &         ! cosine of the half angle of reflection
    !                 perpDistanceOrig, & ! perpendicular distance between neutron path and origin
    !                 sinGamma, &         ! sine of the half angle of reflection
    !                 perpDistance, &     ! distance of a point outside surface to reflective plane
    !                 dOrigin             ! distance from the starting point to the origin

    ! Reflective transforms will not be allowed to occur in geometries other than planes
    call fatalError(Here,'Spheres may not have reflective boundaries')

    ! Calculate the normalised origin vector
    !Ovector= self % origin - r
    !dOrigin=norm2(Ovector)
    !Ovector = Ovector / dOrigin

    ! Obtain a perpendicular vector by taking the cross product
    ! Rotation will take place around this axis
    !perpVector = crossProduct(Ovector,u)

    ! Calculate sinGamma using the cross-product of the two vectors and the sine law
    !sinGamma = (dOrigin/self%radius)*norm2(perpVector)
    !cosGamma = sqrt(1 - sinGamma*sinGamma)

    ! Rotate the particle direction by gamma to obtain the reflection plane normal
    ! Use Rodrigues' rotation formula
    ! (Neglect [1 - cosGamma] term in formula due to taking a dotProduct of perp. vectors)
    !normal = cosGamma*u + sinGamma*crossProduct(perpVector,u)

    ! Obtain the point at which the intersection occurs from the normal, radius and origin
    !intersect = self%origin + normal * self%radius

    ! Calculate particle position outside the sphere
    !r = r + u*distance
    !MAY NEED TO FIX

    ! Calculate the perpendicular distance to the reflecting plane
    !perpDistance = abs(dotProduct(normal,r-intersect))

    ! Translate the particle position across the plane
    !r = r - 2*perpDistance*normal

    ! Reflect the particle direction (independent of intersection point for plane)
    !u = u - 2*dotProduct(normal,u)*normal

  end subroutine reflectiveTransform

  !!
  !! Supply the normal vector given a point on the sphere
  !!
  function normalVector(self, r) result(normal)
    class(sphere), intent(in)               :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3)             :: normal

    normal = TWO * (r - self % origin)

  end function normalVector

  !!
  !! Give an error: this routine should not be called for a non-compound surface
  !!
  function whichSurface(self, r, u) result(surfPointer)
    class(sphere), intent(in)               :: self
    real(defReal), dimension(3), intent(in) :: r, u
    class(surface), pointer                 :: surfPointer
    character(100), parameter :: Here = 'whichSurface (sphere_class.f90)'

    call fatalError(Here,'This function should never be called for a simple surface')

  end function whichSurface

  !!
  !! Set boundary conditions for a sphere: may only be vacuum
  !!
  subroutine setBoundaryConditions(self, BC)
    class(sphere), intent(inout)                :: self
    integer(shortInt), dimension(:), intent(in) :: BC
    character(100), parameter :: Here = 'setBoundaryConditions (sphere_class.f90)'

    if (any(BC /= vacuum)) then
      call fatalError(Here,'Sphere boundaries may only be vacuum')

    else
      self % isVacuum = .TRUE.

    end if
  end subroutine setBoundaryConditions

  !!
  !! Apply boundary conditions
  !!
  subroutine boundaryTransform(self, r, u, isVacuum)
    class(sphere), intent(in)                  :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u
    logical(defBool), intent(inout)            :: isVacuum
    character(100), parameter :: Here = 'boundaryTransform (sphere_class.f90)'

    if (self % isVacuum) then
      isVacuum = .TRUE.

    else
      call fatalError(Here,'This routine should only be called if there sphere has vacuum boundaries')

    end if
  end subroutine boundaryTransform

end module sphere_class
