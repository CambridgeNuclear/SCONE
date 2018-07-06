module xSquareCylinder_class

  use numPrecision
  use universalVariables
  use genericProcedures, only : fatalError
  use dictionary_class,  only : dictionary

  use surface_inter,     only : surface, printSurfDef
  use yPlane_class,      only : yPlane
  use zPlane_class,      only : zPlane

  implicit none
  private

  !!
  !! Constants describing surface properties
  !!
  character(nameLen),parameter :: TYPE_NAME    = 'xSquareCylinder'
  logical(defBool),parameter   :: ACCEPTS_BC   = .true.

  !!
  !! Constructor
  !!
  interface xSquareCylinder
    module procedure xSquareCylinder_fromDict
  end interface

  !!
  !! Square cylinder aligned with x-axis
  !!
  type, public, extends(surface) :: xSquareCylinder
    private
    type(yPlane), dimension(:), pointer :: yPlanes => null()
    type(zPlane), dimension(:), pointer :: zPlanes => null()
    real(defReal), dimension(3)         :: a                ! the half-width in each direction of the cylinder (0 in x)
    real(defReal), dimension(3)         :: origin = [ZERO, ZERO, ZERO]

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

  end type xSquareCylinder


contains

  !!
  !! Initialise the square cylinder from components
  !!
  subroutine init(self, origin, a, id, name)
    class(xSquareCylinder), intent(inout)   :: self
    real(defReal), dimension(3), intent(in) :: a
    real(defReal), dimension(3), intent(in) :: origin
    integer(shortInt), intent(in), optional :: id
    character(*), optional, intent(in)      :: name
    integer(shortInt)                       :: i
    real(defReal)                           :: neg

    self % isCompound = .true.
    self % a = a
    if((a(2) < surface_tol) .OR. (a(3) < surface_tol)) &
    call fatalError('init, xSquareCylinder','Widths must be greater than surface tolerance')
    self % origin = origin

    if(present(id)) self % id = id
    if(present(name)) self % name = name

    if(associated(self % yPlanes)) deallocate (self % yPlanes)
    if(associated(self % zPlanes)) deallocate (self % zPlanes)

    allocate(self%yPlanes(2))
    allocate(self%zPlanes(2))

    ! Initialise each plane in each cardinal direction
    neg = +ONE
    do i = 1, 2
      call self%yPlanes(i)%init(origin(2) + neg*a(2))
      call self%zPlanes(i)%init(origin(3) + neg*a(3))
      neg = -ONE
    end do

  end subroutine init

  !!
  !! Returns and initialised instance of xSquareCylinder from dictionary and name
  !!
  function xSquareCylinder_fromDict(dict,name) result(new)
    class(dictionary), intent(in)  :: dict
    character(nameLen), intent(in) :: name
    type(xSquareCylinder)          :: new
    integer(shortInt)              :: id
    real(defReal), dimension(3)    :: origin, halfwidth
    character(100),parameter :: Here ='xSquareCylinder_fromDict ( xSquareCylinder_class.f90)'

    id = dict % getInt('id')
    if(id < 1) call fatalError(Here,'Invalid surface id provided')

    halfwidth = dict % getRealArray('halfwidth')
    origin = dict % getRealArray('origin')

    call new % init(origin, halfwidth, id, name)

  end function xSquareCylinder_fromDict

  !!
  !! Evaluate the surface function of the square cylinder
  !!
  function evaluate(self, r) result(res)
    class(xSquareCylinder), intent(in)      :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: res
    real(defReal), dimension(2)             :: resY      ! Results from yPlanes
    real(defReal), dimension(2)             :: resZ      ! Results from zPlanes
    real(defReal)                           :: absMinRes ! Identify if particle sits on a surface
    real(defReal)                           :: testRes   ! The most positive residual

    ! Evaluate the front planes' surface functions
    resY(1) = self % yPlanes(1) % evaluate(r)
    resZ(1) = self % zPlanes(1) % evaluate(r)
    absMinRes = min(abs(resY(1)),abs(resZ(1)))
    testRes = max(resY(1),resZ(1))

    ! If the largest result is greater than the surface tolerance
    ! then the particle is not inside the box and the testRes should
    ! be returned
    if (testRes > surface_tol) then
      res = testRes
      return
    end if

    ! Evaluate the rear planes' surface functions
    ! These results are negated to satisfy the definitions
    ! of 'inside' and 'outside'
    resY(2) = -self % yPlanes(2) % evaluate(r)
    resZ(2) = -self % zPlanes(2) % evaluate(r)
    absMinRes = min(absMinRes,abs(resY(2)),abs(resZ(2)))
    testRes = max(resY(2),resZ(2))

    ! If the testRes value is greater than the surface tolerance
    ! then the particle is outside the box
    if (testRes > surface_tol) then
      res = testRes
      return

    ! The particle is either in or on the box
    ! If the absolute minimum value is less than the surface tolerance
    ! then the particle is on the surface
    else if (absMinRes < surface_tol) then
      res = absMinRes
      return

    ! Otherwise, the particle must be inside the box
    ! res will be a negative value
    else
      res = testRes
    end if

  end function evaluate

  !!
  !! Return parameter character containing TYPE NAME
  !!
  function type(self)
    class(xSquareCylinder), intent(in) :: self
    character(nameLen)                 :: type

    type = TYPE_NAME

  end function type

  !!
  !! Return string with definition of this surface
  !!
  pure subroutine getDef(self,string)
    class(xSquareCylinder), intent(in)     :: self
    character(:),allocatable,intent(inout) :: string

    string = printSurfDef(TYPE_NAME, [self % a, self % origin])

  end subroutine getDef

  !!
  !! Override base type function to returns .false.
  !!
  function cannotBeBoundary(self) result(itCant)
    class(xSquareCylinder), intent(in) :: self
    logical(defBool)                   :: itCant

    itCant = .not.ACCEPTS_BC

  end function cannotBeBoundary

  !!
  !! Calculate the distance to the nearest surface of the box
  !! Requires checking that only real surfaces are intercepted,
  !! i.e., not the extensions of the box plane surfaces
  !!
  function distanceToSurface(self, r, u) result(distance)
    class(xSquareCylinder), intent(in)      :: self
    real(defReal), dimension(3), intent(in) :: r, u
    real(defReal), dimension(3)             :: posBound, negBound, testPoint
    real(defReal)                           :: distance
    real(defReal)                           :: testDistance

    distance = INFINITY
    ! Find the positive and negative bounds which the particle
    ! must fall within
    posBound = self % origin + self % a + surface_tol
    negBound = self % origin - self % a - surface_tol

    testDistance = self%yPlanes(1)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    if ((testPoint(3) < posBound(3)) .and. (testPoint(3) > negBound(3))) then
      if (testDistance < distance) distance = testDistance
    end if

    testDistance = self%yPlanes(2)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    if ((testPoint(3) < posBound(3)) .and. (testPoint(3) > negBound(3))) then
      if (testDistance < distance) distance = testDistance
    end if

    testDistance = self%zPlanes(1)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    if ((testPoint(2) < posBound(2)) .and. (testPoint(2) > negBound(2))) then
      if (testDistance < distance) distance = testDistance
    end if

    testDistance = self%zPlanes(2)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    if ((testPoint(2) < posBound(2)) .and. (testPoint(2) > negBound(2))) then
      if (testDistance < distance) distance = testDistance
    end if

  end function distanceToSurface

  !!
  !! Determine on which surface the particle is located and obtain
  !! its normal vector
  !!
  function normalVector(self, r)result(normal)
    class(xSquareCylinder), intent(in)      :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3)             :: normal, posBound, negBound
    character(100),parameter :: Here ='normalVector ( xSquareCylinder_class.f90)'

    ! Compare the point's position to the maximum and minimum
    posBound = self % origin + self % a
    negBound = self % origin - self % a

    if (abs(posBound(2) - r(2)) < surface_tol) then
      normal = self % yPlanes(1) % normalVector(r)
      return
    else if (abs(negBound(2) - r(2)) < surface_tol) then
      normal = self % yPlanes(2) % normalVector(r)
      return
    else if (abs(posBound(3) - r(3)) < surface_tol) then
      normal = self % zPlanes(1) % normalVector(r)
      return
    else if (abs(negBound(3) - r(3)) < surface_tol) then
      normal = self % zPlanes(2) % normalVector(r)
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
    class(xSquareCylinder), intent(in)      :: self
    real(defReal), dimension(3), intent(in) :: r, u
    class(surface), pointer                 :: surfPointer
    real(defReal), dimension(3)             :: testPoint, posBound, negBound
    real(defReal)                           :: testDistance, distance

    distance = INFINITY
    posBound = self % origin + self % a + surface_tol
    negBound = self % origin - self % a - surface_tol

    ! Evaluate distance to each plane and point to surface
    ! with the minimum real distance
    testDistance = self%yPlanes(1)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    if ((testPoint(3) < posBound(3)).and.(testPoint(3) > negBound(3))) then
      if (testDistance < distance) then
        distance = testDistance
        surfPointer => self % yPlanes(1)
      end if
    end if

    testDistance = self%yPlanes(2)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    if ((testPoint(3) < posBound(3)) .and. (testPoint(3) > negBound(3))) then
      if (testDistance < distance) then
        distance = testDistance
        surfPointer => self % yPlanes(2)
      end if
    end if

    testDistance = self%zPlanes(1)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    if ((testPoint(2) < posBound(2)) .and. (testPoint(2) > negBound(2))) then
      if (testDistance < distance) then
        distance = testDistance
        surfPointer => self % zPlanes(1)
      end if
    end if

    testDistance = self%zPlanes(2)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    if ((testPoint(2) < posBound(2)) .and. (testPoint(2) > negBound(2))) then
      if (testDistance < distance) then
        distance = testDistance
        surfPointer => self % zPlanes(2)
      end if
    end if

  end function whichSurface

  !!
  !! Set boundary conditions for an xSquareCylinder: may only be vacuum
  !!
  subroutine setBoundaryConditions(self, BC)
    class(xSquareCylinder), intent(inout)       :: self
    integer(shortInt), dimension(:), intent(in) :: BC
    character(100),parameter :: Here ='setBoundaryConditions ( xSquareCylinder_class.f90)'

    if(size(BC) < 6) call fatalError(Here,'Wrong size of BC string. Must be at least 6')

    ! Positive y boundary
    if(BC(3) == vacuum) then
      self % yPlanes(1) % isVacuum = .TRUE.
    else if(BC(3) == reflective) then
      self % yPlanes(1) % isReflective = .TRUE.
    else if(BC(3) == periodic) then
      if(BC(4) /= periodic) then
        call fatalError(Here, 'Both positive and negative boundary conditions must be periodic')

      else
        self % yPlanes(1) % isPeriodic = .TRUE.
        self % yPlanes(1) % periodicTranslation = [ZERO, -TWO*self % a(2), ZERO]
      end if
    else
      call fatalError(Here,'Invalid boundary condition provided')
    end if

    ! Negative y boundary
    if(BC(4) == vacuum) then
      self % yPlanes(2) % isVacuum = .TRUE.
    else if(BC(4) == reflective) then
      self % yPlanes(2) % isReflective = .TRUE.
    else if(BC(4) == periodic) then
      if(BC(3) /= periodic) then
        call fatalError(Here, 'Both positive and negative boundary conditions must be periodic')
      else
        self % yPlanes(2) % isPeriodic = .TRUE.
        self % yPlanes(2) % periodicTranslation = [ZERO, TWO*self % a(2), ZERO]
      end if
    else
      call fatalError(Here,'Invalid boundary condition provided')

    end if

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
        self % zPlanes(1) % periodicTranslation = [ZERO, ZERO, -2.*self % a(3)]
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
        self % zPlanes(2) % periodicTranslation = [ZERO, ZERO, 2.*self % a(3)]
      end if
    else
      call fatalError(Here,'Invalid boundary condition provided')
    end if
  end subroutine setBoundaryConditions

  !!
  !! Apply boundary transformation
  !!
  subroutine boundaryTransform(self, r, u, isVacuum)
    class(xSquareCylinder), intent(in)         :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u
    logical(defBool), intent(inout)            :: isVacuum
    logical(defBool)                           :: left, right, above, below
    character(100),parameter :: Here ='boundaryTransform ( xSquareCylinder_class.f90)'

    ! Point can be within one of 9 regions
    ! Locate which region and then apply BCs as appropriate
    left = self % yPlanes(1) % halfspace(r, u)
    right = .NOT. self % yPlanes(2) % halfspace(r, u)
    above = self % zPlanes(1) % halfspace(r, u)
    below = .NOT. self % zPlanes(2) % halfspace(r, u)

    if (left .AND. above) then
      call self % yPlanes(1) % boundaryTransform(r, u, isVacuum)
      call self % zPlanes(1) % boundaryTransform(r, u, isVacuum)
    else if (right .AND. above) then
      call self % yPlanes(2) % boundaryTransform(r, u, isVacuum)
      call self % zPlanes(1) % boundaryTransform(r, u, isVacuum)
    else if (above) then
      call self % zPlanes(1) % boundaryTransform(r, u, isVacuum)
    else if (left .AND. below) then
      call self % yPlanes(1) % boundaryTransform(r, u, isVacuum)
      call self % zPlanes(2) % boundaryTransform(r, u, isVacuum)
    else if (right .AND. below) then
      call self % yPlanes(2) % boundaryTransform(r, u, isVacuum)
      call self % zPlanes(2) % boundaryTransform(r, u, isVacuum)
    else if (below) then
      call self % zPlanes(2) % boundaryTransform(r, u, isVacuum)
    else if (left) then
      call self % yPlanes(1) % boundaryTransform(r, u, isVacuum)
    else if (right) then
      call self % yPlanes(2) % boundaryTransform(r, u, isVacuum)
    else
      call fatalError(Here,'Cannot apply boundary condition: point is inside surface')
    end if

  end subroutine boundaryTransform

end module xSquareCylinder_class
