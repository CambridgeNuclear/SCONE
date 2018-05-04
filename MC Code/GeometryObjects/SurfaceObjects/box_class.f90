!
! Derived surface type: cuboid composed of 6 plane surfaces
!
module box_class
  use numPrecision
  use genericProcedures, only : fatalError
  use universalVariables

  use surface_class
  use plane_class

  implicit none
  private

  type, public, extends(surface) :: box
    type(xPlane), dimension(:), pointer :: xPlanes => null()
    type(yPlane), dimension(:), pointer :: yPlanes => null()
    type(zPlane), dimension(:), pointer :: zPlanes => null()
    real(defReal), dimension(3)         :: a    ! the half-width in each direction of the cubic box
    real(defReal), dimension(3)         :: origin = [ZERO, ZERO, ZERO]
  contains
    procedure :: init => initBox
    procedure :: evaluate => evaluateBox
    procedure :: distanceToSurface => distanceToBox
    procedure :: reflectiveTransform => reflectiveTransformBox
    procedure :: normalVector => normalVectorBox
    procedure :: whichSurface => whichPlane
    procedure :: setBoundaryConditions => setBoundaryConditionsBox
    procedure :: boundaryTransform => boundaryTransformBox
  end type box

contains

  !
  ! Initialise the box as six plane surfaces
  ! Employ the convention that the first plane
  ! is in front of the origin, the second behind
  !
  subroutine initBox(self, origin, a, id, name)
    class(box), intent(inout)               :: self
    real(defReal), dimension(3), intent(in) :: a, &
                                              origin
    integer(shortInt), intent(in), optional :: id
    character(*), optional, intent(in)      :: name
    integer(shortInt)                       :: i
    real(defReal)                           :: neg

    self % isCompound = .true.
    self % origin = origin
    self % a = a
    if(any(a < surface_tol)) &
    call fatalError('initBox, box','Box dimensions must be greater than surface tolerance')
    if(present(id)) self % id = id
    if (present(name)) self % name = name

    if(associated(self % xPlanes)) deallocate (self % xPlanes)
    if(associated(self % yPlanes)) deallocate (self % yPlanes)
    if(associated(self % zPlanes)) deallocate (self % zPlanes)

    allocate(self % xPlanes(2))
    allocate(self % yPlanes(2))
    allocate(self % zPlanes(2))

    ! Initialise each plane in each cardinal direction
    neg = +ONE
    do i = 1, 2
      call self % xPlanes(i) % init(origin(1) + neg*a(1))
      call self % yPlanes(i) % init(origin(2) + neg*a(2))
      call self % zPlanes(i) % init(origin(3) + neg*a(3))
      neg = -ONE
    end do

  end subroutine initBox

  !
  ! Evaluate the surface function of the box
  !
  function evaluateBox(self, r) result(res)
    class(box), intent(in)                  :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: res
    real(defReal), dimension(2)             :: resX      ! Results from xPlanes
    real(defReal), dimension(2)             :: resY      ! Results from yPlanes
    real(defReal), dimension(2)             :: resZ      ! Results from zPlanes
    real(defReal)                           :: absMinRes ! Identify if particle sits on a surface
    real(defReal)                           :: testRes   ! The most positive residual

    ! Evaluate the front planes' surface functions
    resX(1) = self % xPlanes(1) % evaluate(r)
    resY(1) = self % yPlanes(1) % evaluate(r)
    resZ(1) = self % zPlanes(1) % evaluate(r)
    absMinRes = min(abs(resX(1)),abs(resY(1)),abs(resZ(1)))
    testRes = max(resX(1),resY(1),resZ(1))

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
    resZ(2) = -self % zPlanes(2) % evaluate(r)
    absMinRes = min(absMinRes,abs(resX(2)),abs(resY(2)),abs(resZ(2)))
    testRes = max(resX(2),resY(2),resZ(2))

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

  end function evaluateBox

  !
  ! Calculate the distance to the nearest surface of the box
  ! Requires checking that only real surfaces are intercepted,
  ! i.e., not the extensions of the box plane surfaces
  !
  function distanceToBox(self, r, u) result(distance)
    class(box), intent(in)                  :: self
    real(defReal), dimension(3), intent(in) :: r, u
    real(defReal), dimension(3)             :: posBound, negBound, testPoint
    real(defReal)                           :: distance
    real(defReal)                           :: testDistance

    distance = INFINITY
    ! Find the positive and negative bounds which the particle
    ! must fall within
    posBound = self % origin + self % a
    negBound = self % origin - self % a

    testDistance = self%xPlanes(1)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    if ((testPoint(2) < posBound(2)) .and. (testPoint(2) > negBound(2)) .and. &
    (testPoint(3) < posBound(3)) .and. (testPoint(3) > negBound(3))) then
      if (testDistance < distance) distance = testDistance
    end if

    testDistance = self%xPlanes(2)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    if ((testPoint(2) < posBound(2)) .and. (testPoint(2) > negBound(2)) .and. &
    (testPoint(3) < posBound(3)) .and. (testPoint(3) > negBound(3))) then
      if (testDistance < distance) distance = testDistance
    end if

    testDistance = self%yPlanes(1)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    if ((testPoint(1) < posBound(1)) .and. (testPoint(1) > negBound(1)) .and. &
    (testPoint(3) < posBound(3)) .and. (testPoint(3) > negBound(3))) then
      if (testDistance < distance) distance = testDistance
    end if

    testDistance = self%yPlanes(2)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    if ((testPoint(1) < posBound(1)) .and. (testPoint(1) > negBound(1)) .and. &
    (testPoint(3) < posBound(3)) .and. (testPoint(3) > negBound(3))) then
      if (testDistance < distance) distance = testDistance
    end if

    testDistance = self%zPlanes(1)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    if ((testPoint(2) < posBound(2)) .and. (testPoint(2) > negBound(2)) .and. &
    (testPoint(1) < posBound(1)) .and. (testPoint(1) > negBound(1))) then
      if (testDistance < distance) distance = testDistance
    end if

    testDistance = self%zPlanes(2)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    if ((testPoint(2) < posBound(2)) .and. (testPoint(2) > negBound(2)) .and. &
    (testPoint(1) < posBound(1)) .and. (testPoint(1) > negBound(1))) then
      if (testDistance < distance) distance = testDistance
    end if

  end function distanceToBox

  !
  ! Apply a reflective transformation to a particle during delta tracking
  ! Assumes particle has already crossed the plane, hence the direction reversal
  ! Do so by determining which plane the particle intersects and applying the plane reflection
  !
  ! This routine has been obviated due to BC implementation in the transport operator and surface
  ! search in the cell
  !
  subroutine reflectiveTransformBox(self, r, u)
    class(box), intent(in)                     :: self
    real(defReal), dimension(3), intent(inout) :: r, u
    class(surface), pointer                    :: surfPointer

    surfPointer => self % whichSurface(r, u)
    call surfPointer % reflectiveTransform(r,u)

  end subroutine reflectiveTransformBox

  !
  ! Determine on which surface the particle is located and obtain
  ! its normal vector
  !
  function normalVectorBox(self, r) result(normal)
    class(box), intent(in)                  :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3)             :: normal
    real(defReal), dimension(3)             :: posBound, negBound

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
    else if (abs(posBound(3) - r(3)) < surface_tol) then
      normal = self % zPlanes(1) % normalVector(r)
      return
    else if (abs(negBound(3) - r(3)) < surface_tol) then
      normal = self % zPlanes(2) % normalVector(r)
      return
    else
      call fatalError('normalVector, box','Point is not on a surface')
    end if

  end function normalVectorBox

  !
  ! Helper routine: find the plane which a particle will have crossed
  ! This is done by calculating the distance to each surface
  ! Include inside or outside to allow generality (need to change
  ! particle direction otherwise)
  !
  function whichPlane(self, r, u) result(surfPointer)
    class(box), intent(in)                  :: self
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
    if ((testPoint(2) < posBound(2)).and.(testPoint(2) > negBound(2)) .and. &
    (testPoint(3) < posBound(3)).and.(testPoint(3) > negBound(3))) then
      if (testDistance < distance) then
        distance = testDistance
        surfPointer => self % xPlanes(1)
      end if
    end if

    testDistance = self%xPlanes(2)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    if ((testPoint(2) < posBound(2)).and.(testPoint(2) > negBound(2)) .and. &
    (testPoint(3) < posBound(3)).and.(testPoint(3) > negBound(3))) then
      if (testDistance < distance) then
        distance = testDistance
        surfPointer => self % xPlanes(2)
      end if
    end if

    testDistance = self%yPlanes(1)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    if ((testPoint(1) < posBound(1)).and.(testPoint(1) > negBound(1)) .and. &
    (testPoint(3) < posBound(3)).and.(testPoint(3) > negBound(3))) then
      if (testDistance < distance) then
        distance = testDistance
        surfPointer => self % yPlanes(1)
      end if
    end if

    testDistance = self%yPlanes(2)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    if ((testPoint(1) < posBound(1)) .and. (testPoint(1) > negBound(1)) .and. &
    (testPoint(3) < posBound(3)) .and. (testPoint(3) > negBound(3))) then
      if (testDistance < distance) then
        distance = testDistance
        surfPointer => self % yPlanes(2)
      end if
    end if

    testDistance = self%zPlanes(1)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    if ((testPoint(2) < posBound(2)) .and. (testPoint(2) > negBound(2)) .and. &
    (testPoint(1) < posBound(1)) .and. (testPoint(1) > negBound(1))) then
      if (testDistance < distance) then
        distance = testDistance
        surfPointer => self % zPlanes(1)
      end if
    end if

    testDistance = self%zPlanes(2)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    if ((testPoint(2) < posBound(2)) .and. (testPoint(2) > negBound(2)) .and. &
    (testPoint(1) < posBound(1)) .and. (testPoint(1) > negBound(1))) then
      if (testDistance < distance) then
        distance = testDistance
        surfPointer => self % zPlanes(2)
      end if
    end if

  end function whichPlane

  !!
  !! Apply generic boundary conditions to the box
  !!
  subroutine setBoundaryConditionsBox(self, BC)
    class(box), intent(inout)                   :: self
    integer(shortInt), dimension(6), intent(in) :: BC

    ! Positive x boundary
    if(BC(1) == vacuum) then
      self % xPlanes(1) % isVacuum = .TRUE.
    else if(BC(1) == reflective) then
      self % xPlanes(1) % isReflective = .TRUE.
    else if(BC(1) == periodic) then
      if(BC(2) /= periodic) then
        call fatalError('setBoundaryConditions, box', &
        'Both positive and negative boundary conditions must be periodic')
      else
        self % xPlanes(1) % isPeriodic = .TRUE.
        self % xPlanes(1) % periodicTranslation = [-TWO*self % a(1), ZERO, ZERO]
      end if
    else
      call fatalError('setBoundaryConditionsBox','Invalid boundary condition provided')
    end if

    ! Negative x boundary
    if(BC(2) == vacuum) then
      self % xPlanes(2) % isVacuum = .TRUE.
    else if(BC(2) == reflective) then
      self % xPlanes(2) % isReflective = .TRUE.
    else if(BC(2) == periodic) then
      if(BC(1) /= periodic) then
        call fatalError('setBoundaryConditions, box', &
        'Both positive and negative boundary conditions must be periodic')
      else
        self % xPlanes(2) % isPeriodic = .TRUE.
        self % xPlanes(2) % periodicTranslation = [TWO*self % a(1), ZERO, ZERO]
      end if
    else
      call fatalError('setBoundaryConditionsBox','Invalid boundary condition provided')
    end if

    ! Positive y boundary
    if(BC(3) == vacuum) then
      self % yPlanes(1) % isVacuum = .TRUE.
    else if(BC(3) == reflective) then
      self % yPlanes(1) % isReflective = .TRUE.
    else if(BC(3) == periodic) then
      if(BC(4) /= periodic) then
        call fatalError('setBoundaryConditions, box', &
        'Both positive and negative boundary conditions must be periodic')
      else
        self % yPlanes(1) % isPeriodic = .TRUE.
        self % yPlanes(1) % periodicTranslation = [ZERO, -TWO*self % a(2), ZERO]
      end if
    else
      call fatalError('setBoundaryConditionsBox','Invalid boundary condition provided')
    end if

    ! Negative y boundary
    if(BC(4) == vacuum) then
      self % yPlanes(2) % isVacuum = .TRUE.
    else if(BC(4) == reflective) then
      self % yPlanes(2) % isReflective = .TRUE.
    else if(BC(4) == periodic) then
      if(BC(3) /= periodic) then
        call fatalError('setBoundaryConditions, box', &
        'Both positive and negative boundary conditions must be periodic')
      else
        self % yPlanes(2) % isPeriodic = .TRUE.
        self % yPlanes(2) % periodicTranslation = [ZERO, TWO*self % a(2), ZERO]
      end if
    else
      call fatalError('setBoundaryConditionsBox','Invalid boundary condition provided')
    end if

    ! Positive z boundary
    if(BC(5) == vacuum) then
      self % zPlanes(1) % isVacuum = .TRUE.
    else if(BC(5) == reflective) then
      self % zPlanes(1) % isReflective = .TRUE.
    else if(BC(5) == periodic) then
      if(BC(6) /= periodic) then
        call fatalError('setBoundaryConditions, box', &
        'Both positive and negative boundary conditions must be periodic')
      else
        self % zPlanes(1) % isPeriodic = .TRUE.
        self % zPlanes(1) % periodicTranslation = [ZERO, ZERO, -TWO*self % a(3)]
      end if
    else
      call fatalError('setBoundaryConditionsBox','Invalid boundary condition provided')
    end if

    ! Negative z boundary
    if(BC(6) == vacuum) then
      self % zPlanes(2) % isVacuum = .TRUE.
    else if(BC(6) == reflective) then
      self % zPlanes(2) % isReflective = .TRUE.
    else if(BC(6) == periodic) then
      if(BC(5) /= periodic) then
        call fatalError('setBoundaryConditions, box', &
        'Both positive and negative boundary conditions must be periodic')
      else
        self % zPlanes(2) % isPeriodic = .TRUE.
        self % zPlanes(2) % periodicTranslation = [ZERO, ZERO, TWO*self % a(3)]
      end if
    else
      call fatalError('setBoundaryConditionsBox','Invalid boundary condition provided')
    end if

  end subroutine setBoundaryConditionsBox

  !!
  !! Check halfspaces to identify which combination of
  !! boundary conditions to apply to a point
  !!
  subroutine boundaryTransformBox(self, r, u, isVacuum)
    class(box), intent(in)                     :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u
    logical(defBool), intent(inout)            :: isVacuum
    logical(defBool)                           :: front, back, left, right, above, below

    ! Point can be within one of 27 regions
    ! Locate which region and then apply BCs as appropriate
    front = self % xPlanes(1) % halfspace(r, u)
    back = .NOT. self % xPlanes(2) % halfspace(r, u)
    left = self % yPlanes(1) % halfspace(r, u)
    right = .NOT. self % yPlanes(2) % halfspace(r, u)
    above = self % zPlanes(1) % halfspace(r, u)
    below = .NOT. self % zPlanes(2) % halfspace(r, u)

    ! Point is in the upper front left corner
    if (front .AND. left .AND. above) then
      call self % xPlanes(1) % boundaryTransform(r, u, isVacuum)
      call self % yPlanes(1) % boundaryTransform(r, u, isVacuum)
      call self % zPlanes(1) % boundaryTransform(r, u, isVacuum)

    ! Point is in the upper back left corner
    else if (back .AND. left .AND. above) then
      call self % xPlanes(2) % boundaryTransform(r, u, isVacuum)
      call self % yPlanes(1) % boundaryTransform(r, u, isVacuum)
      call self % zPlanes(1) % boundaryTransform(r, u, isVacuum)

    ! Point is above and to the left
    else if (left .AND. above) then
      call self % yPlanes(1) % boundaryTransform(r, u, isVacuum)
      call self % zPlanes(1) % boundaryTransform(r, u, isVacuum)

    ! Point is in the upper front right corner
    else if (front .AND. right .AND. above) then
      call self % xPlanes(1) % boundaryTransform(r, u, isVacuum)
      call self % yPlanes(2) % boundaryTransform(r, u, isVacuum)
      call self % zPlanes(1) % boundaryTransform(r, u, isVacuum)

    ! Point is in the upper back right corner
    else if (back .AND. right .AND. above) then
      call self % xPlanes(2) % boundaryTransform(r, u, isVacuum)
      call self % yPlanes(2) % boundaryTransform(r, u, isVacuum)
      call self % zPlanes(1) % boundaryTransform(r, u, isVacuum)

    ! Point is above and to the right
    else if (right .AND. above) then
      call self % yPlanes(2) % boundaryTransform(r, u, isVacuum)
      call self % zPlanes(1) % boundaryTransform(r, u, isVacuum)

    ! Point is above
    else if (above) then
      call self % zPlanes(1) % boundaryTransform(r, u, isVacuum)

    ! Point is above and in front
    else if (front .AND. above) then
      call self % xPlanes(1) % boundaryTransform(r, u, isVacuum)
      call self % zPlanes(1) % boundaryTransform(r, u, isVacuum)

    ! Point is above and behind
    else if (back .AND. above) then
      call self % xPlanes(2) % boundaryTransform(r, u, isVacuum)
      call self % zPlanes(1) % boundaryTransform(r, u, isVacuum)

    ! Point is in the lower front left
    else if (front .AND. left .AND. below) then
      call self % xPlanes(1) % boundaryTransform(r, u, isVacuum)
      call self % yPlanes(1) % boundaryTransform(r, u, isVacuum)
      call self % zPlanes(2) % boundaryTransform(r, u, isVacuum)

    ! Point is in the lower back left
    else if (back .AND. left .AND. below) then
      call self % xPlanes(2) % boundaryTransform(r, u, isVacuum)
      call self % yPlanes(1) % boundaryTransform(r, u, isVacuum)
      call self % zPlanes(2) % boundaryTransform(r, u, isVacuum)

    ! Point is in the lower left
    else if (left .AND. below) then
      call self % yPlanes(1) % boundaryTransform(r, u, isVacuum)
      call self % zPlanes(2) % boundaryTransform(r, u, isVacuum)

    ! Point is in the lower front right
    else if (front .AND. right .AND. below) then
      call self % xPlanes(1) % boundaryTransform(r, u, isVacuum)
      call self % yPlanes(2) % boundaryTransform(r, u, isVacuum)
      call self % zPlanes(2) % boundaryTransform(r, u, isVacuum)

    ! Point is in the lower back right
    else if (back .AND. right .AND. below) then
      call self % xPlanes(2) % boundaryTransform(r, u, isVacuum)
      call self % yPlanes(2) % boundaryTransform(r, u, isVacuum)
      call self % zPlanes(2) % boundaryTransform(r, u, isVacuum)

    ! Point is in the lower right
    else if (right .AND. below) then
      call self % yPlanes(2) % boundaryTransform(r, u, isVacuum)
      call self % zPlanes(2) % boundaryTransform(r, u, isVacuum)

    ! Point is below
    else if (below) then
      call self % zPlanes(2) % boundaryTransform(r, u, isVacuum)

    ! Point is in the front left
    else if (front .AND. left) then
      call self % xPlanes(1) % boundaryTransform(r, u, isVacuum)
      call self % yPlanes(1) % boundaryTransform(r, u, isVacuum)

    ! Point is in the back left
    else if (back .AND. left) then
      call self % xPlanes(2) % boundaryTransform(r, u, isVacuum)
      call self % yPlanes(1) % boundaryTransform(r, u, isVacuum)

    ! Point is on the left
    else if (left) then
      call self % yPlanes(1) % boundaryTransform(r, u, isVacuum)

    ! Point is on the front right
    else if (front .AND. right) then
      call self % xPlanes(1) % boundaryTransform(r, u, isVacuum)
      call self % yPlanes(2) % boundaryTransform(r, u, isVacuum)

    ! Point is on the back right
    else if (back .AND. right) then
      call self % xPlanes(2) % boundaryTransform(r, u, isVacuum)
      call self % yPlanes(2) % boundaryTransform(r, u, isVacuum)

    ! Point is on the front
    else if (front) then
      call self % xPlanes(1) % boundaryTransform(r, u, isVacuum)

    ! Point is on the back
    else if (back) then
      call self % xPlanes(2) % boundaryTransform(r, u, isVacuum)

    ! Point is in the box
    else
      call fatalError('boundaryTransformation, box',&
      'Cannot apply boundary transformation: point is already within the boundary')
    end if

  end subroutine boundaryTransformBox

end module box_class
