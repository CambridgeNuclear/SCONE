module ySquareCylinder_class

  use numPrecision
  use genericProcedures, only : fatalError
  use universalVariables

  use surface_inter,     only : surface
  use xPlane_class,      only : xPlane
  use zPlane_class,      only : zPlane

  implicit none
  private

  type, public, extends(surface) :: ySquareCylinder
    private
    type(xPlane), dimension(:), pointer :: xPlanes => null()
    type(zPlane), dimension(:), pointer :: zPlanes => null()
    real(defReal), dimension(3)         :: a                ! the half-width in each direction of the cylinder (0 in y)
    real(defReal), dimension(3)         :: origin = [ZERO, ZERO, ZERO]

  contains
    procedure :: init
    procedure :: evaluate
    procedure :: distanceToSurface
    procedure :: normalVector
    procedure :: reflectiveTransform
    procedure :: whichSurface
    procedure :: setBoundaryConditions
    procedure :: boundaryTransform

  end type ySquareCylinder

contains
  !!
  !! Initialise the cylinder as four plane surfaces
  !!
  subroutine init(self, origin, a, id, name)
    class(ySquareCylinder), intent(inout)   :: self
    real(defReal), dimension(3), intent(in) :: a
    real(defReal), dimension(3), intent(in) :: origin
    integer(shortInt), intent(in), optional :: id
    character(*), optional, intent(in)      :: name
    integer(shortInt)                       :: i
    real(defReal)                           :: neg

    self % isCompound = .true.
    self % a = a
    if((a(1) < surface_tol) .OR. (a(3) < surface_tol)) &
    call fatalError('init, ySquareCylinder','Widths must be greater than surface tolerance')
    self % origin = origin

    if(present(id)) self % id = id
    if(present(name)) self % name = name

    if(associated(self % xPlanes)) deallocate (self % xPlanes)
    if(associated(self % zPlanes)) deallocate (self % zPlanes)

    allocate(self%xPlanes(2))
    allocate(self%zPlanes(2))

    ! Initialise each plane in each cardinal direction
    neg = +ONE
    do i = 1, 2
      call self%xPlanes(i)%init(origin(1) + neg*a(1))
      call self%zPlanes(i)%init(origin(3) + neg*a(3))
      neg = -ONE
    end do

  end subroutine init

  !!
  !! Evaluate the surface function of the square cylinder
  !!
  function evaluate(self, r) result(res)
    class(ySquareCylinder), intent(in)      :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: res
    real(defReal), dimension(2)             :: resX      ! Results from yPlanes
    real(defReal), dimension(2)             :: resZ      ! Results from zPlanes
    real(defReal)                           :: absMinRes ! Identify if particle sits on a surface
    real(defReal)                           :: testRes   ! The most positive residual

    ! Evaluate the front planes' surface functions
    resX(1) = self % xPlanes(1) % evaluate(r)
    resZ(1) = self % zPlanes(1) % evaluate(r)
    absMinRes = min(abs(resX(1)),abs(resZ(1)))
    testRes = max(resX(1),resZ(1))

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
    resZ(2) = -self % zPlanes(2) % evaluate(r)
    absMinRes = min(absMinRes,abs(resX(2)),abs(resZ(2)))
    testRes = max(resX(2),resZ(2))

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
  !! Calculate the distance to the nearest surface of the cylinder
  !! Requires checking that only real surfaces are intercepted,
  !! i.e., not the extensions of the plane surfaces
  !!
  function distanceToSurface(self, r, u) result(distance)
    class(ySquareCylinder), intent(in)      :: self
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
    if ((testPoint(3) < posBound(3)) .and. (testPoint(3) > negBound(3))) then
      if (testDistance < distance) distance = testDistance
    end if

    testDistance = self%xPlanes(2)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    if ((testPoint(3) < posBound(3)) .and. (testPoint(3) > negBound(3))) then
      if (testDistance < distance) distance = testDistance
    end if

    testDistance = self%zPlanes(1)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    if ((testPoint(1) < posBound(1)) .and. (testPoint(1) > negBound(1))) then
      if (testDistance < distance) distance = testDistance
    end if

    testDistance = self%zPlanes(2)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    if ((testPoint(1) < posBound(1)) .and. (testPoint(1) > negBound(1))) then
      if (testDistance < distance) distance = testDistance
    end if

  end function distanceToSurface

  !!
  !! Apply a reflective transformation to a particle during delta tracking
  !! Do so by determining which plane the particle intersects and applying the plane reflection
  !!
  !! This routine has been obviated due to BC implementation in the transport operator and surface
  !! search in the cell
  !!
  subroutine reflectiveTransform(self, r, u)
    class(ySquareCylinder), intent(in)         :: self
    real(defReal), dimension(3), intent(inout) :: r, u
    class(surface), pointer                    :: surfPointer

    call fatalError('reflectiveTransformYSquCyl','This routine should not be called')
    surfPointer => self % whichSurface(r, u)
    call surfPointer % reflectiveTransform(r,u)

  end subroutine reflectiveTransform

  !!
  !! Determine on which surface the particle is located and obtain
  !! its normal vector
  !!
  function normalVector(self, r) result(normal)
    class(ySquareCylinder), intent(in)      :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3)             :: normal, posBound, negBound

    ! Compare the point's position to the maximum and minimum
    posBound = self % origin + self % a
    negBound = self % origin - self % a

    if (abs(posBound(1) - r(1)) < surface_tol) then
      normal = self % xPlanes(1) % normalVector(r)
      return
    else if (abs(negBound(1) - r(1)) < surface_tol) then
      normal = self % xPlanes(2) % normalVector(r)
      return
    else if (abs(posBound(3) - r(3)) < surface_tol) then
      normal = self % zPlanes(1) % normalVector(r)
      return
    else if (abs(negBound(3) - r(3)) < surface_tol) then
      normal = self % zPlanes(2) % normalVector(r)
      return
    else
      call fatalError('normalVector, ySquareCylinder','Point is not on a surface')
    end if

  end function normalVector

  !!
  !! Helper routine: find the plane which a particle will have crossed
  !! This is done by calculating the distance to each surface
  !! Include inside or outside to allow generality (need to change
  !! particle direction otherwise)
  !!
  function whichSurface(self, r, u) result(surfPointer)
    class(ySquareCylinder), intent(in)      :: self
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
    if ((testPoint(3) < posBound(3)).and.(testPoint(3) > negBound(3))) then
      if (testDistance < distance) then
        distance = testDistance
        surfPointer => self % xPlanes(1)
      end if
    end if

    testDistance = self%xPlanes(2)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    if ((testPoint(3) < posBound(3)) .and. (testPoint(3) > negBound(3))) then
      if (testDistance < distance) then
        distance = testDistance
        surfPointer => self % xPlanes(2)
      end if
    end if

    testDistance = self%zPlanes(1)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    if ((testPoint(1) < posBound(1)) .and. (testPoint(1) > negBound(1))) then
      if (testDistance < distance) then
        distance = testDistance
        surfPointer => self % zPlanes(1)
      end if
    end if

    testDistance = self%zPlanes(2)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    if ((testPoint(1) < posBound(1)) .and. (testPoint(1) > negBound(1))) then
      if (testDistance < distance) then
        distance = testDistance
        surfPointer => self % zPlanes(2)
      end if
    end if

  end function whichSurface

  !!
  !! Apply generic boundary conditions to a ySquareCylinder
  !!
  subroutine setBoundaryConditions(self, BC)
    class(ySquareCylinder), intent(inout)       :: self
    integer(shortInt), dimension(6), intent(in) :: BC

    ! Positive x boundary
    if(BC(1) == vacuum) then
      self % xPlanes(1) % isVacuum = .TRUE.
    else if(BC(1) == reflective) then
      self % xPlanes(1) % isReflective = .TRUE.
    else if(BC(1) == periodic) then
      if(BC(2) /= periodic) then
        call fatalError('setBoundaryConditionsYSquareCylinder', &
        'Both positive and negative boundary conditions must be periodic')
      else
        self % xPlanes(1) % isPeriodic = .TRUE.
        self % xPlanes(1) % periodicTranslation = [-TWO*self % a(1), ZERO, ZERO]
      end if
    else
      call fatalError('setBoundaryConditionsYSquareCylinder','Invalid boundary condition provided')
    end if

    ! Negative x boundary
    if(BC(2) == vacuum) then
      self % xPlanes(2) % isVacuum = .TRUE.
    else if(BC(2) == reflective) then
      self % xPlanes(2) % isReflective = .TRUE.
    else if(BC(2) == periodic) then
      if(BC(1) /= periodic) then
        call fatalError('setBoundaryConditionsYSquareCylinder', &
        'Both positive and negative boundary conditions must be periodic')
      else
        self % xPlanes(2) % isPeriodic = .TRUE.
        self % xPlanes(2) % periodicTranslation = [TWO*self % a(1), ZERO, ZERO]
      end if
    else
      call fatalError('setBoundaryConditionsYSquareCylinder','Invalid boundary condition provided')
    end if

    ! Positive z boundary
    if(BC(5) == vacuum) then
      self % zPlanes(1) % isVacuum = .TRUE.
    else if(BC(5) == reflective) then
      self % zPlanes(1) % isReflective = .TRUE.
    else if(BC(5) == periodic) then
      if(BC(6) /= periodic) then
        call fatalError('setBoundaryConditionsYSquareCylinder', &
        'Both positive and negative boundary conditions must be periodic')
      else
        self % zPlanes(1) % isPeriodic = .TRUE.
        self % zPlanes(1) % periodicTranslation = [ZERO, ZERO, -TWO*self % a(3)]
      end if
    else
      call fatalError('setBoundaryConditionsYSquareCylinder','Invalid boundary condition provided')
    end if

    ! Negative z boundary
    if(BC(6) == vacuum) then
      self % zPlanes(2) % isVacuum = .TRUE.
    else if(BC(6) == reflective) then
      self % zPlanes(2) % isReflective = .TRUE.
    else if(BC(6) == periodic) then
      if(BC(5) /= periodic) then
        call fatalError('setBoundaryConditionsYSquareCylinder', &
        'Both positive and negative boundary conditions must be periodic')
      else
        self % zPlanes(2) % isPeriodic = .TRUE.
        self % zPlanes(2) % periodicTranslation = [ZERO, ZERO, TWO*self % a(3)]
      end if
    else
      call fatalError('setBoundaryConditionsYSquareCylinder','Invalid boundary condition provided')
    end if

  end subroutine setBoundaryConditions

  !!
  !! Apply boundary transformation
  !!
  subroutine boundaryTransform(self, r, u, isVacuum)
    class(ySquareCylinder), intent(in)         :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u
    logical(defBool), intent(inout)            :: isVacuum
    logical(defBool)                           :: front, back, above, below

    ! Point can be within one of 9 regions
    ! Locate which region and then apply BCs as appropriate
    front = self % xPlanes(1) % halfspace(r, u)
    back = .NOT. self % xPlanes(2) % halfspace(r, u)
    above = self % zPlanes(1) % halfspace(r, u)
    below = .NOT. self % zPlanes(2) % halfspace(r, u)

    if (front .AND. above) then
      call self % zPlanes(1) % boundaryTransform(r, u, isVacuum)
      call self % zPlanes(1) % boundaryTransform(r, u, isVacuum)
    else if (back .AND. above) then
      call self % xPlanes(2) % boundaryTransform(r, u, isVacuum)
      call self % zPlanes(1) % boundaryTransform(r, u, isVacuum)
    else if (above) then
      call self % zPlanes(1) % boundaryTransform(r, u, isVacuum)
    else if (front .AND. below) then
      call self % xPlanes(1) % boundaryTransform(r, u, isVacuum)
      call self % zPlanes(2) % boundaryTransform(r, u, isVacuum)
    else if (back .AND. below) then
      call self % xPlanes(2) % boundaryTransform(r, u, isVacuum)
      call self % zPlanes(2) % boundaryTransform(r, u, isVacuum)
    else if (below) then
      call self % zPlanes(2) % boundaryTransform(r, u, isVacuum)
    else if (front) then
      call self % xPlanes(1) % boundaryTransform(r, u, isVacuum)
    else if (back) then
      call self % xPlanes(2) % boundaryTransform(r, u, isVacuum)
    else
      call fatalError('boundaryTransform, ySquareCylinder',&
      'Cannot apply boundary condition: point is inside surface')
    end if

  end subroutine boundaryTransform

end module ySquareCylinder_class
