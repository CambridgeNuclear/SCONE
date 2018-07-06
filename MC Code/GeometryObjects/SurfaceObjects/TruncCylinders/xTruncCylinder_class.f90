module xTruncCylinder_class

  use numPrecision
  use universalVariables
  use genericProcedures, only : fatalError
  use dictionary_class,  only : dictionary

  use surface_inter,     only : surface, printSurfDef
  use xPlane_class,      only : xPlane
  use xCylinder_class,   only : xCylinder

  implicit none
  private

  !!
  !! Constants describing surface properties
  !!
  character(nameLen),parameter :: TYPE_NAME    = 'xTruncCylinder'
  logical(defBool),parameter   :: ACCEPTS_BC   = .true.

  !!
  !! Constructor
  !!
  interface xTruncCylinder
    module procedure xTruncCylinder_fromDict
  end interface

  !!
  !! Truncated cylinder aligned with x-axis
  !! Has finate length
  !!
  type, public, extends(surface) :: xTruncCylinder
    private
    type(xPlane), dimension(:), pointer   :: xPlanes => null()
    type(xCylinder), pointer              :: cyl => null()  ! substituent cylinder
    real(defReal)                         :: a              ! the half-width separation of the planes
    real(defReal), dimension(3)           :: origin = [ZERO, ZERO, ZERO]
  contains
    procedure :: init
    procedure :: evaluate
    procedure :: type
    procedure :: getDef
    procedure :: cannotBeBoundary
    procedure :: distanceToSurface
    procedure :: normalVector
    procedure :: whichSurface
    procedure :: setBoundaryConditions
    procedure :: boundaryTransform

  end type xTruncCylinder

contains

  !!
  !! Initialise the box as six plane surfaces
  !!
  subroutine init(self, origin, a, radius, id, name)
    class(xTruncCylinder), intent(inout)     :: self
    real(defReal), dimension(3), intent(in)  :: origin
    real(defReal), intent(in)                :: a
    real(defReal), intent(in)                :: radius
    integer(shortInt), intent(in), optional  :: id
    character(*), optional, intent(in)       :: name

    self % isCompound = .true.
    self % a = a
    if(a < surface_tol) &
    call fatalError('init, xTruncCylinder','Width must be greater than surface tolerance')
    if(radius < surface_tol) &
    call fatalError('init, xTruncCylinder','Radius must be greater than surface tolerance')
    self % origin = origin
    if(present(id)) self % id = id
    if(present(name)) self % name = name

    if(associated(self % xPlanes)) deallocate (self % xPlanes)

    allocate(self%xPlanes(2))

    call self % xPlanes(1) % init(origin(1) + a)
    call self % xPlanes(2) % init(origin(1) - a)
    call self % cyl % init(radius, origin)

    ! Hash and store surface definition
    call self % hashSurfDef()

  end subroutine init

  !!
  !! Returns and initialised instance of xTruncCylinder from dictionary and name
  !!
  function xTruncCylinder_fromDict(dict,name) result(new)
    class(dictionary), intent(in)  :: dict
    character(nameLen), intent(in) :: name
    type(xTruncCylinder)           :: new
    integer(shortInt)              :: id
    real(defReal)                  :: radius, halfheight
    real(defReal), dimension(3)    :: origin
    character(100),parameter :: Here ='xTruncCylinder_fromDict ( xTruncCylinder_class.f90)'

    id = dict % getInt('id')
    if(id < 1) call fatalError(Here,'Invalid surface id provided')

    radius = dict % getReal('radius')
    halfheight = dict % getReal('halfheight')
    origin = dict % getRealArray('origin')

    call new % init(origin, halfheight, radius, id, name)

  end function xTruncCylinder_fromDict

  !!
  !! Evaluate the surface function of the square cylinder
  !!
  function evaluate(self, r) result(res)
    class(xTruncCylinder), intent(in)       :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: res
    real(defReal), dimension(2)             :: resX      ! Results from yPlanes
    real(defReal)                           :: resCyl    ! Results from cylinder
    real(defReal)                           :: absMinRes ! Identify if particle sits on a surface
    real(defReal)                           :: testRes   ! The most positive residual

    ! Evaluate the top plane's surface function
    resX(1) = self % xPlanes(1) % evaluate(r)

    ! If the result is greater than the surface tolerance
    ! then the particle is not inside the cylinder and the result should
    ! be returned
    if (resX(1) > surface_tol) then
      res = resX(1)
      return
    end if

    resX(2) = -self % xPlanes(2) % evaluate(r)
    absMinRes = min(abs(resX(1)),abs(resX(2)))
    testRes = max(resX(1),resX(2))

    ! Same as before - particle is outside the cylinder
    if (resX(2) > surface_tol) then
      res = resX(2)
      return
    end if

    ! Evaluate the cylinder's surface function
    resCyl = self % cyl % evaluate(r)
    absMinRes = min(absMinRes,abs(resCyl))
    testRes = max(testRes, resCyl)

    ! If the cylinder result is greater than the surface tolerance
    ! then the particle is outside the cylinder
    if (resCyl > surface_tol) then
      res = resCyl
      return

    ! The particle is either in or on the cylinder
    ! If the absolute minimum value is less than the surface tolerance
    ! then the particle is on the surface
    else if (absMinRes < surface_tol) then
      res = absMinRes
      return

    ! The particle is inside the cylinder
    ! res should be a negative value
    else
      res = testRes
      return
    end if

  end function evaluate

  !!
  !! Return parameter character containing TYPE NAME
  !!
  function type(self)
    class(xTruncCylinder), intent(in) :: self
    character(nameLen)     :: type

    type = TYPE_NAME

  end function type

  !!
  !! Return string with definition of this surface
  !!
  pure subroutine getDef(self,string)
    class(xTruncCylinder), intent(in)      :: self
    character(:),allocatable,intent(inout) :: string

    ! This chained call is SUPER dirty...
    string = printSurfDef(TYPE_NAME, [self % cyl % radius, self % a, self % origin])

  end subroutine getDef

  !!
  !! Override base type function to returns .false.
  !!
  function cannotBeBoundary(self) result(itCant)
    class(xTruncCylinder), intent(in) :: self
    logical(defBool)       :: itCant

    itCant = .not.ACCEPTS_BC

  end function cannotBeBoundary

  !!
  !! Calculate the distance to the nearest surface of the cylinder
  !! Requires checking that only real surfaces are intercepted,
  !! i.e., not the extensions of the plane or cylinderical surfaces
  !!
  function distanceToSurface(self, r, u)result(distance)
    class(xTruncCylinder), intent(in)       :: self
    real(defReal), dimension(3), intent(in) :: r, u
    real(defReal), dimension(3)             :: testPoint
    real(defReal)                           :: posBound, negBound
    real(defReal)                           :: distance
    real(defReal)                           :: testDistance
    real(defReal)                           :: testCyl

    distance = INFINITY
    ! Find the positive and negative bounds which the particle
    ! must fall within
    posBound = self % xPlanes(1) % x0
    negBound = self % xPlanes(2) % x0

    testDistance = self%xPlanes(1)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    testCyl = self % cyl % evaluate(testPoint)

    ! Ensure point is within the radius of the cylinder
    if (testCyl < surface_tol)  then
      if (testDistance < distance) distance = testDistance
    end if

    testDistance = self%xPlanes(2)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    testCyl = self % cyl % evaluate(testPoint)

    if (testCyl < surface_tol) then
      if (testDistance < distance) distance = testDistance
    end if

    testDistance = self% cyl %distanceToSurface(r,u)
    testPoint = r + u*testDistance
    if ((testPoint(1) < posBound) .and. (testPoint(1) > negBound)) then
      if (testDistance < distance) distance = testDistance
    end if

  end function distanceToSurface

  !!
  !! Determine on which surface the particle is located and obtain
  !! its normal vector
  !!
  function normalVector(self, r) result(normal)
    class(xTruncCylinder), intent(in)       :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3)             :: normal
    real(defReal)                           :: posBound, negBound
    character(100),parameter :: Here ='normalVector ( xTruncCylinder_class.f90)'

    ! Compare the point's position to the maximum and minimum
    posBound = self % xPlanes(1) % x0
    negBound = self % xPlanes(2) % x0

    if (abs(posBound - r(1)) < surface_tol) then
      normal = self % xPlanes(1) % normalVector(r)
      return
    else if (abs(negBound - r(1)) < surface_tol) then
      normal = self % xPlanes(2) % normalVector(r)
      return
    else if (abs(self % cyl % evaluate(r)) < surface_tol) then
      normal = self % cyl % normalVector(r)
      return
    else
      call fatalError(Here,'Point is not on a surface')
    end if

  end function normalVector

  !!
  !! Helper routine: find the plane which a particle will have crossed
  !! This is done by calculating the distance to each surface
  !! Include inside or outside to allow generality (need to change
  !! particle direction otherwise)
  !!
  function whichSurface(self, r, u) result(surfPointer)
    class(xTruncCylinder), intent(in)       :: self
    real(defReal), dimension(3), intent(in) :: r, u
    class(surface), pointer                 :: surfPointer
    real(defReal), dimension(3)             :: testPoint
    real(defReal)                           :: posBound, negBound
    real(defReal)                           :: testDistance, distance
    real(defReal)                           :: testCyl

    distance = INFINITY
    posBound = self % xPlanes(1) % x0
    negBound = self % xPlanes(2) % x0

    ! Evaluate distance to each plane and point to surface
    ! with the minimum real distance
    testDistance = self%xPlanes(1)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    testCyl = self % cyl % evaluate(testPoint)

    if (testCyl < surface_tol) then
      if (testDistance < distance) then
        distance = testDistance
        surfPointer => self % xPlanes(1)
      end if
    end if

    testDistance = self%xPlanes(2)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    testCyl = self % cyl % evaluate(testPoint)

    if (testCyl < surface_tol) then
      if (testDistance < distance) then
        distance = testDistance
        surfPointer => self % xPlanes(2)
      end if
    end if

    testDistance = self % cyl % distanceToSurface(r,u)
    testPoint = r + u*testDistance
    if ((testPoint(1) < posBound).and.(testPoint(1) > negBound))  then
      if (testDistance < distance) then
        distance = testDistance
        surfPointer => self % cyl
      end if
    end if

  end function whichSurface

  !!
  !! Set boundary conditions for an xTruncCylinder
  !!
  subroutine setBoundaryConditions(self, BC)
    class(xTruncCylinder), intent(inout)        :: self
    integer(shortInt), dimension(:), intent(in) :: BC
    character(100),parameter :: Here ='setBoundaryConditions ( xTruncCylinder_class.f90)'

    if(size(BC) < 6) call fatalError(Here,'Wrong size of BC string. Must be at least 6')

    ! Positive x boundary
    if(BC(1) == vacuum) then
      self % xPlanes(1) % isVacuum = .TRUE.
    else if(BC(1) == reflective) then
      self % xPlanes(1) % isReflective = .TRUE.
    else if(BC(1) == periodic) then
      if(BC(2) /= periodic) then
        call fatalError(Here, 'Both positive and negative boundary conditions must be periodic')
      else
        self % xPlanes(1) % isPeriodic = .TRUE.
        self % xPlanes(1) % periodicTranslation = [-TWO*self % a, ZERO, ZERO]
      end if
    else
      call fatalError(Here,'Invalid boundary condition provided')
    end if

    ! Negative x boundary
    if(BC(2) == vacuum) then
      self % xPlanes(2) % isVacuum = .TRUE.
    else if(BC(2) == reflective) then
      self % xPlanes(2) % isReflective = .TRUE.
    else if(BC(2) == periodic) then
      if(BC(1) /= periodic) then
        call fatalError(Here, 'Both positive and negative boundary conditions must be periodic')
      else
        self % xPlanes(2) % isPeriodic = .TRUE.
        self % xPlanes(2) % periodicTranslation = [TWO*self % a, ZERO, ZERO]
      end if
    else
      call fatalError(Here,'Invalid boundary condition provided')
    end if

    if(any(BC(3:6) /= vacuum))then
      call fatalError(Here,'Cylinder boundaries may only be vacuum')
    else
      self % cyl % isVacuum = .TRUE.
    end if
  end subroutine setBoundaryConditions

  !!
  !! Apply boundary transformation
  !!
  subroutine boundaryTransform(self, r, u, isVacuum)
    class(xTruncCylinder), intent(in)          :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u
    logical(defBool), intent(inout)            :: isVacuum
    logical(defBool)                           :: front, back, outsideCyl
    character(100),parameter :: Here ='boundaryTransform ( xTruncCylinder_class.f90)'

    outsideCyl = self % cyl % halfspace(r, u)
    if (outsideCyl) then
      call self % cyl % boundaryTransform(r, u, isVacuum)
      return
    end if

    front = self % xPlanes(1) % halfspace(r, u)
    back = .NOT. self % xPlanes(2) % halfspace(r, u)
    if (front) then
      call self % xPlanes(1) % boundaryTransform(r, u, isVacuum)
    else if (back) then
      call self % xPlanes(2) % boundaryTransform(r, u, isVacuum)
    else
      call fatalError(Here,'Cannot apply boundary condition: point is inside the surface')
    end if

  end subroutine boundaryTransform

end module xTruncCylinder_class
