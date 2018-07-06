module zSquareCylinder_class

  use numPrecision
  use universalVariables
  use genericProcedures, only : fatalError
  use dictionary_class,  only : dictionary

  use surface_inter,     only : surface, printSurfDef
  use xPlane_class,      only : xPlane
  use yPlane_class,      only : yPlane

  implicit none
  private

  !!
  !! Constants describing surface properties
  !!
  character(nameLen),parameter :: TYPE_NAME    = 'zSquareCylinder'
  logical(defBool),parameter   :: ACCEPTS_BC   = .true.

  !!
  !! Constructor
  !!
  interface zSquareCylinder
    module procedure zSquareCylinder_fromDict
  end interface

  !!
  !! Square cylinder aligned with y-axis
  !!
  type, public, extends(surface) :: zSquareCylinder
    private
    type(xPlane), dimension(:), pointer :: xPlanes => null()
    type(yPlane), dimension(:), pointer :: yPlanes => null()
    real(defReal), dimension(3)         :: a                 ! the half-width in each direction of the cylinder (0 in z)
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

  end type zSquareCylinder


contains

  !!
  !! Initialise the square cylinder from components
  !!
  subroutine init(self, origin, a, id, name)
    class(zSquareCylinder), intent(inout)   :: self
    real(defReal), dimension(3), intent(in) :: a
    real(defReal), dimension(3), intent(in) :: origin
    integer(shortInt), intent(in), optional :: id
    character(*), optional, intent(in)      :: name
    integer(shortInt)                       :: i
    real(defReal)                           :: neg

    self % isCompound = .true.
    self % a = a
    if((a(1) < surface_tol) .OR. (a(2) < surface_tol)) &
    call fatalError('init, zSquareCylinder','Widths must be greater than surface tolerance')
    self % origin = origin

    if(present(id)) self % id = id
    if(present(name)) self % name = name

    if(associated(self % xPlanes)) deallocate (self % xPlanes)
    if(associated(self % yPlanes)) deallocate (self % yPlanes)

    allocate(self%xPlanes(2))
    allocate(self%yPlanes(2))

    ! Initialise each plane in each cardinal direction
    neg = +ONE
    do i = 1, 2
      call self%xPlanes(i)%init(origin(1) + neg*a(1))
      call self%yPlanes(i)%init(origin(2) + neg*a(2))
      neg = -ONE
    end do

    ! Hash and store surface definition
    call self % hashSurfDef()

  end subroutine init

  !!
  !! Returns and initialised instance of ySquareCylinder from dictionary and name
  !!
  function zSquareCylinder_fromDict(dict,name) result(new)
    class(dictionary), intent(in)  :: dict
    character(nameLen), intent(in) :: name
    type(zSquareCylinder)          :: new
    integer(shortInt)              :: id
    real(defReal), dimension(3)    :: origin, halfwidth
    character(100),parameter :: Here ='zSquareCylinder_fromDict ( zSquareCylinder_class.f90)'

    id = dict % getInt('id')
    if(id < 1) call fatalError(Here,'Invalid surface id provided')

    halfwidth = dict % getRealArray('halfwidth')
    origin = dict % getRealArray('origin')

    call new % init(origin, halfwidth, id, name)

  end function zSquareCylinder_fromDict

  !!
  !! Evaluate the surface function of the square cylinder
  !!
  function evaluate(self, r) result(res)
    class(zSquareCylinder), intent(in)      :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: res
    real(defReal), dimension(2)             :: resX      ! Results from yPlanes
    real(defReal), dimension(2)             :: resY      ! Results from zPlanes
    real(defReal)                           :: absMinRes ! Identify if particle sits on a surface
    real(defReal)                           :: testRes   ! The most positive residual

    ! Evaluate the front planes' surface functions
    resX(1) = self % xPlanes(1) % evaluate(r)
    resY(1) = self % yPlanes(1) % evaluate(r)
    absMinRes = min(abs(resX(1)),abs(resY(1)))
    testRes = max(resX(1),resY(1))

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
    resX(2) = -self % xPlanes(2) % evaluate(r)
    resY(2) = -self % yPlanes(2) % evaluate(r)
    absMinRes = min(absMinRes,abs(resX(2)),abs(resY(2)))
    testRes = max(resX(2),resY(2))

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
    class(zSquareCylinder), intent(in) :: self
    character(nameLen)                 :: type

    type = TYPE_NAME

  end function type

  !!
  !! Return string with definition of this surface
  !!
  pure subroutine getDef(self,string)
    class(zSquareCylinder), intent(in)     :: self
    character(:),allocatable,intent(inout) :: string

    string = printSurfDef(TYPE_NAME, [self % a, self % origin])

  end subroutine getDef

  !!
  !! Override base type function to returns .false.
  !!
  function cannotBeBoundary(self) result(itCant)
    class(zSquareCylinder), intent(in) :: self
    logical(defBool)                   :: itCant

    itCant = .not.ACCEPTS_BC

  end function cannotBeBoundary

  !!
  !! Calculate the distance to the nearest surface of the cylinder
  !! Requires checking that only real surfaces are intercepted,
  !! i.e., not the extensions of the plane surfaces
  !!
  function distanceToSurface(self, r, u)result(distance)
    class(zSquareCylinder), intent(in)      :: self
    real(defReal), dimension(3), intent(in) :: r, u
    real(defReal), dimension(3)             :: posBound, negBound, testPoint
    real(defReal)                           :: distance
    real(defReal)                           :: testDistance

    distance = INFINITY
    ! Find the positive and negative bounds which the particle
    ! must fall within
    posBound = self % origin + self % a + surface_tol
    negBound = self % origin - self % a - surface_tol

    testDistance = self%xPlanes(1)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    if ((testPoint(2) < posBound(2)) .and. (testPoint(2) > negBound(2))) then
      if (testDistance < distance) distance = testDistance
    end if

    testDistance = self%xPlanes(2)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    if ((testPoint(2) < posBound(2)) .and. (testPoint(2) > negBound(2))) then
      if (testDistance < distance) distance = testDistance
    end if

    testDistance = self%yPlanes(1)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    if ((testPoint(1) < posBound(1)) .and. (testPoint(1)> negBound(1))) then
      if (testDistance < distance) distance = testDistance
    end if

    testDistance = self%yPlanes(2)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    if ((testPoint(1) < posBound(1)) .and. (testPoint(1) > negBound(1))) then
      if (testDistance < distance) distance = testDistance
    end if

  end function distanceToSurface

  !!
  !! Determine on which surface the particle is located and obtain
  !! its normal vector
  !!
  function normalVector(self, r) result(normal)
    class(zSquareCylinder), intent(in)      :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3)             :: normal, posBound, negBound
    character(100),parameter :: Here ='normalVector ( zSquareCylinder_class.f90)'

    ! Compare the point's position to the maximum and minimum
    posBound = self % origin + self % a
    negBound = self % origin - self % a

    if (abs(posBound(1) - r(1)) < surface_tol) then
      normal = self % xPlanes(1) % normalVector(r)
      return
    else if (abs(negBound(1) - r(1)) < surface_tol) then
      normal = self % xPlanes(2) % normalVector(r)
      return
    else if (abs(posBound(2) - r(2)) < surface_tol) then
      normal = self % yPlanes(1) % normalVector(r)
      return
    else if (abs(negBound(2) - r(2)) < surface_tol) then
      normal = self % yPlanes(2) % normalVector(r)
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
    class(zSquareCylinder), intent(in)      :: self
    real(defReal), dimension(3), intent(in) :: r, u
    class(surface), pointer                 :: surfPointer
    real(defReal), dimension(3)             :: testPoint, posBound, negBound
    real(defReal)                           :: testDistance, distance

    distance = INFINITY
    posBound = self % origin + self % a + surface_tol
    negBound = self % origin - self % a - surface_tol

    ! Evaluate distance to each plane and point to surface
    ! with the minimum real distance
    testDistance = self%xPlanes(1)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    if (((testPoint(2) < posBound(2)).and.(testPoint(2) > negBound(2)))) then
      if (testDistance < distance) then
        distance = testDistance
        surfPointer => self % xPlanes(1)
      end if
    end if

    testDistance = self%xPlanes(2)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    if (((testPoint(2) < posBound(2)) .and. (testPoint(2) > negBound(2)))) then
      if (testDistance < distance) then
        distance = testDistance
        surfPointer => self % xPlanes(2)
      end if
    end if

    testDistance = self%yPlanes(1)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    if (((testPoint(1) < posBound(1)) .and. (testPoint(1) > negBound(1)))) then
      if (testDistance < distance) then
        distance = testDistance
        surfPointer => self % yPlanes(1)
      end if
    end if

    testDistance = self%yPlanes(2)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    if (((testPoint(1) < posBound(1)) .and. (testPoint(1) > negBound(1)))) then
      if (testDistance < distance) then
        distance = testDistance
        surfPointer => self % yPlanes(2)
      end if
    end if

  end function whichSurface

  !!
  !! Apply generic boundary conditions to a zSquareCylinder
  !!
  subroutine setBoundaryConditions(self, BC)
    class(zSquareCylinder), intent(inout)       :: self
    integer(shortInt), dimension(:), intent(in) :: BC
    character(100),parameter :: Here ='setBoundaryConditions ( zSquareCylinder_class.f90)'

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
        self % xPlanes(1) % periodicTranslation = [-TWO*self % a(1), ZERO, ZERO]
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
        self % xPlanes(2) % periodicTranslation = [TWO*self % a(1), ZERO, ZERO]
      end if
    else
      call fatalError(Here,'Invalid boundary condition provided')
    end if

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

  end subroutine setBoundaryConditions

  !!
  !! Apply boundary transformation
  !!
  subroutine boundaryTransform(self, r, u, isVacuum)
    class(zSquareCylinder), intent(in)         :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u
    logical(defBool), intent(inout)            :: isVacuum
    logical(defBool)                           :: left, right, front, back
    character(100),parameter :: Here ='boundaryTransform ( zSquareCylinder_class.f90)'

    ! Point can be within one of 9 regions
    ! Locate which region and then apply BCs as appropriate
    left = self % yPlanes(1) % halfspace(r, u)
    right = .NOT. self % yPlanes(2) % halfspace(r, u)
    front = self % xPlanes(1) % halfspace(r, u)
    back = .NOT. self % xPlanes(2) % halfspace(r, u)

    if (left .AND. front) then
      call self % yPlanes(1) % boundaryTransform(r, u, isVacuum)
      call self % xPlanes(1) % boundaryTransform(r, u, isVacuum)
    else if (right .AND. front) then
      call self % yPlanes(2) % boundaryTransform(r, u, isVacuum)
      call self % xPlanes(1) % boundaryTransform(r, u, isVacuum)
    else if (front) then
      call self % xPlanes(1) % boundaryTransform(r, u, isVacuum)
    else if (left .AND. back) then
      call self % yPlanes(1) % boundaryTransform(r, u, isVacuum)
      call self % xPlanes(2) % boundaryTransform(r, u, isVacuum)
    else if (right .AND. back) then
      call self % yPlanes(2) % boundaryTransform(r, u, isVacuum)
      call self % xPlanes(2) % boundaryTransform(r, u, isVacuum)
    else if (back) then
      call self % xPlanes(2) % boundaryTransform(r, u, isVacuum)
    else if (left) then
      call self % yPlanes(1) % boundaryTransform(r, u, isVacuum)
    else if (right) then
      call self % yPlanes(2) % boundaryTransform(r, u, isVacuum)
    else
      call fatalError(Here,'Cannot apply boundary condition: point is inside surface')
    end if

  end subroutine boundaryTransform

end module zSquareCylinder_class
