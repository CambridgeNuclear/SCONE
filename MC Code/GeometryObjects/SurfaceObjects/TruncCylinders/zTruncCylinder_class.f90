module zTruncCylinder_class

  use numPrecision
  use universalVariables
  use genericProcedures, only : fatalError
  use dictionary_class,  only : dictionary

  use surface_inter,     only : surface, printSurfDef
  use zPlane_class,      only : zPlane
  use zCylinder_class,   only : zCylinder

  implicit none
  private

  !!
  !! Constants describing surface properties
  !!
  character(nameLen),parameter :: TYPE_NAME    = 'zTruncCylinder'
  logical(defBool),parameter   :: ACCEPTS_BC   = .true.

  !!
  !! Constructor
  !!
  interface zTruncCylinder
    module procedure zTruncCylinder_fromDict
  end interface

  !!
  !! Truncated cylinder aligned with z-axis
  !! Has finate length
  !!
  type, public, extends(surface) :: zTruncCylinder
    private
    type(zPlane), dimension(:), pointer   :: zPlanes => null()
    type(zCylinder), pointer              :: cyl => null()  ! substituent cylinder
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

  end type zTruncCylinder

contains

  !!
  !! Initialise the box as six plane surfaces
  !!
  subroutine init(self, origin, a, radius, id, name)
    class(zTruncCylinder), intent(inout)     :: self
    real(defReal), dimension(3), intent(in)  :: origin
    real(defReal), intent(in)                :: a, radius
    integer(shortInt), intent(in), optional  :: id
    character(*), optional, intent(in)       :: name

    self % isCompound = .true.
    self % a = a
    if(a < surface_tol) &
    call fatalError('init, zTruncCylinder','Width must be greater than surface tolerance')
    if(radius < surface_tol) &
    call fatalError('init, zTruncCylinder','Radius must be greater than surface tolerance')
    self % origin = origin
    if (present(id)) self % id = id
    if (present(name)) self % name = name

    if(associated(self % zPlanes)) deallocate (self % zPlanes)

    allocate(self%zPlanes(2))

    call self % zPlanes(1) % init(origin(3) + a)
    call self % zPlanes(2) % init(origin(3) - a)
    call self % cyl % init(radius, origin)

    ! Hash and store surface definition
    call self % hashSurfDef()

  end subroutine init

  !!
  !! Returns and initialised instance of zTruncCylinder from dictionary and name
  !!
  function zTruncCylinder_fromDict(dict,name) result(new)
    class(dictionary), intent(in)  :: dict
    character(nameLen), intent(in) :: name
    type(zTruncCylinder)           :: new
    integer(shortInt)              :: id
    real(defReal)                  :: radius, halfheight
    real(defReal), dimension(3)    :: origin
    character(100),parameter :: Here ='zTruncCylinder_fromDict ( zTruncCylinder_class.f90)'

    id = dict % getInt('id')
    if(id < 1) call fatalError(Here,'Invalid surface id provided')

    radius = dict % getReal('radius')
    halfheight = dict % getReal('halfheight')
    origin = dict % getRealArray('origin')

    call new % init(origin, halfheight, radius, id, name)

  end function zTruncCylinder_fromDict


  !!
  !! Evaluate the surface function of the square cylinder
  !!
  function evaluate(self, r) result(res)
    class(zTruncCylinder), intent(in)       :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: res
    real(defReal), dimension(2)             :: resZ     ! Results from zPlanes
    real(defReal)                           :: resCyl   ! Results from cylinder
    real(defReal)                           :: absMinRes! Identify if particle sits on a surface
    real(defReal)                           :: testRes  ! The most positive residual

    ! Evaluate the top plane's surface function
    resZ(1) = self % zPlanes(1) % evaluate(r)

    ! If the result is greater than the surface tolerance
    ! then the particle is not inside the cylinder and the result should
    ! be returned
    if (resZ(1) > surface_tol) then
      res = resZ(1)
      return
    end if

    resZ(2) = -self % zPlanes(2) % evaluate(r)
    absMinRes = min(abs(resZ(1)),abs(resZ(2)))
    testRes = max(resZ(1),resZ(2))

    ! Same as before - particle is outside the cylinder
    if (resZ(2) > surface_tol) then
      res = resZ(2)
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
    class(zTruncCylinder), intent(in) :: self
    character(nameLen)                :: type

    type = TYPE_NAME

  end function type

  !!
  !! Return string with definition of this surface
  !!
  pure subroutine getDef(self,string)
    class(zTruncCylinder), intent(in)      :: self
    character(:),allocatable,intent(inout) :: string

    ! This chained call is SUPER dirty...
    string = printSurfDef(TYPE_NAME, [self % cyl % radius, self % a, self % origin])

  end subroutine getDef

  !!
  !! Override base type function to returns .false.
  !!
  function cannotBeBoundary(self) result(itCant)
    class(zTruncCylinder), intent(in) :: self
    logical(defBool)                  :: itCant

    itCant = .not.ACCEPTS_BC

  end function cannotBeBoundary

  !!
  !! Calculate the distance to the nearest surface of the cylinder
  !! Requires checking that only real surfaces are intercepted,
  !! i.e., not the extensions of the plane or cylinderical surfaces
  !!
  function distanceToSurface(self, r, u) result(distance)
    class(zTruncCylinder), intent(in)       :: self
    real(defReal), dimension(3), intent(in) :: r, u
    real(defReal), dimension(3)             :: testPoint
    real(defReal)                           :: posBound, negBound
    real(defReal)                           :: distance
    real(defReal)                           :: testDistance
    real(defReal)                           :: testCyl

    distance = INFINITY
    ! Find the positive and negative bounds which the particle
    ! must fall within
    posBound = self % zPlanes(1) % z0
    negBound = self % zPlanes(2) % z0

    testDistance = self%zPlanes(1)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    testCyl = self % cyl % evaluate(testPoint)

    ! Ensure point is within the radius of the cylinder
    if (testCyl < surface_tol)  then
      if (testDistance < distance) distance = testDistance
    end if

    testDistance = self%zPlanes(2)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    testCyl = self % cyl % evaluate(testPoint)

    if (testCyl < surface_tol) then
      if (testDistance < distance) distance = testDistance
    end if

    testDistance = self% cyl %distanceToSurface(r,u)
    testPoint = r + u*testDistance
    if ((testPoint(3) < posBound) .and. (testPoint(3) > negBound)) then
      if (testDistance < distance) distance = testDistance
    end if

  end function distanceToSurface

  !!
  !! Determine on which surface the particle is located and obtain
  !! its normal vector
  !!
  function normalVector(self, r) result(normal)
    class(zTruncCylinder), intent(in)       :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3)             :: normal
    real(defReal)                           :: posBound, negBound
    character(100),parameter :: Here ='normalVector ( zTruncCylinder_class.f90)'

    ! Compare the point's position to the maximum and minimum
    posBound = self % zPlanes(1) % z0
    negBound = self % zPlanes(2) % z0

    if (abs(posBound - r(3)) < surface_tol) then
      normal = self % zPlanes(1) % normalVector(r)
      return
    else if (abs(negBound - r(3)) < surface_tol) then
      normal = self % zPlanes(2) % normalVector(r)
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
    class(zTruncCylinder), intent(in)       :: self
    real(defReal), dimension(3), intent(in) :: r, u
    class(surface), pointer                 :: surfPointer
    real(defReal), dimension(3)             :: testPoint
    real(defReal)                           :: posBound, negBound
    real(defReal)                           :: testDistance, distance
    real(defReal)                           :: testCyl

    distance = INFINITY
    posBound = self % zPlanes(1) % z0
    negBound = self % zPlanes(2) % z0

    ! Evaluate distance to each plane and point to surface
    ! with the minimum real distance
    testDistance = self%zPlanes(1)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    testCyl = self % cyl % evaluate(testPoint)

    if (testCyl < surface_tol) then
      if (testDistance < distance) then
        distance = testDistance
        surfPointer => self % zPlanes(1)
      end if
    end if

    testDistance = self%zPlanes(2)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    testCyl = self % cyl % evaluate(testPoint)

    if (testCyl < surface_tol) then
      if (testDistance < distance) then
        distance = testDistance
        surfPointer => self % zPlanes(2)
      end if
    end if

    testDistance = self % cyl % distanceToSurface(r,u)
    testPoint = r + u*testDistance
    if ((testPoint(3) < posBound).and.(testPoint(3) > negBound))  then
      if (testDistance < distance) then
        distance = testDistance
        surfPointer => self % cyl
      end if
    end if

  end function whichSurface

  !!
  !! Set boundary conditions for a zTruncCylinder
  !!
  subroutine setBoundaryConditions(self, BC)
    class(zTruncCylinder), intent(inout)        :: self
    integer(shortInt), dimension(:), intent(in) :: BC
    character(100),parameter :: Here ='setBoundaryConditions ( zTruncCylinder_class.f90)'

    if(size(BC) < 6) call fatalError(Here,'Wrong size of BC string. Must be at least 6')

    ! Positive z boundary
    if(BC(5) == vacuum) then
      self % zPlanes(1) % isVacuum = .TRUE.
    else if(BC(5) == reflective) then
      self % zPlanes(1) % isReflective = .TRUE.
    else if(BC(5) == periodic) then
      if(BC(6) /= periodic) then
        call fatalError(Here, 'Both positive and negative boundary conditions must be periodic')
      else
        self % zPlanes(1) % isPeriodic = .TRUE.
        self % zPlanes(1) % periodicTranslation = [ZERO, ZERO, -TWO*self % a]
      end if
    else
      call fatalError(Here,'Invalid boundary condition provided')
    end if

    ! Negative z boundary
    if(BC(6) == vacuum) then
      self % zPlanes(2) % isVacuum = .TRUE.
    else if(BC(6) == reflective) then
      self % zPlanes(2) % isReflective = .TRUE.
    else if(BC(6) == periodic) then
      if(BC(5) /= periodic) then
        call fatalError(Here, 'Both positive and negative boundary conditions must be periodic')
      else
        self % zPlanes(2) % isPeriodic = .TRUE.
        self % zPlanes(2) % periodicTranslation = [ZERO, ZERO, TWO*self % a]
      end if
    else
      call fatalError(Here,'Invalid boundary condition provided')
    end if

    if(any(BC(1:4) /= vacuum))then
      call fatalError(Here,'Cylinder boundaries may only be vacuum')
    else
      self % cyl % isVacuum = .TRUE.
    end if
  end subroutine setBoundaryConditions

  !!
  !! Apply boundary transformation
  !!
  subroutine boundaryTransform(self, r, u, isVacuum)
    class(zTruncCylinder), intent(in)          :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u
    logical(defBool), intent(inout)            :: isVacuum
    logical(defBool)                           :: above, below, outsideCyl
    character(100),parameter :: Here ='boundaryTransform ( zTruncCylinder_class.f90)'

    outsideCyl = self % cyl % halfspace(r, u)
    if (outsideCyl) then
      call self % cyl % boundaryTransform(r, u, isVacuum)
      return
    end if

    above = self % zPlanes(1) % halfspace(r, u)
    below = .NOT. self % zPlanes(2) % halfspace(r, u)
    if (above) then
      call self % zPlanes(1) % boundaryTransform(r, u, isVacuum)
    else if (below) then
      call self % zPlanes(2) % boundaryTransform(r, u, isVacuum)
    else
      call fatalError(Here,'Cannot apply boundary condition: point is inside the surface')
    end if

  end subroutine boundaryTransform

end module zTruncCylinder_class
