module truncatedCylinder_class
  use numPrecision
  use genericProcedures, only : fatalError
  use universalVariables

  use surface_inter
  use plane_class
  use cylinder_class

  implicit none
  private

  type, public, extends(surface) :: xTruncatedCylinder
    private
    type(xPlane), dimension(:), pointer   :: xPlanes => null()
    type(xCylinder), pointer              :: cyl => null()  ! substituent cylinder
    real(defReal)                         :: a              ! the half-width separation of the planes
    real(defReal), dimension(3)           :: origin = [ZERO, ZERO, ZERO]
  contains
    procedure :: init => initXTruncCyl
    procedure :: evaluate => evaluateXTruncCyl
    procedure :: distanceToSurface => distanceToXTruncCyl
    procedure :: normalVector => normalXTruncCyl
    procedure :: reflectiveTransform => reflectiveTransformXTruncCyl
    procedure :: whichSurface => whichSurfaceXTruncCyl
    procedure :: setBoundaryConditions => setBoundaryConditionsXTruncCylinder
    procedure :: boundaryTransform => boundaryTransformXTruncCylinder
  end type xTruncatedCylinder

  type, public, extends(surface) :: yTruncatedCylinder
    private
    type(yPlane), dimension(:), pointer   :: yPlanes => null()
    type(yCylinder), pointer              :: cyl => null()  ! substituent cylinder
    real(defReal)                         :: a              ! the half-width separation of the planes
    real(defReal), dimension(3)           :: origin = [ZERO, ZERO, ZERO]
  contains
    procedure :: init => initYTruncCyl
    procedure :: evaluate => evaluateYTruncCyl
    procedure :: distanceToSurface => distanceToYTruncCyl
    procedure :: normalVector => normalYTruncCyl
    procedure :: reflectiveTransform => reflectiveTransformYTruncCyl
    procedure :: whichSurface => whichSurfaceYTruncCyl
    procedure :: setBoundaryConditions => setBoundaryConditionsYTruncCylinder
    procedure :: boundaryTransform => boundaryTransformYTruncCylinder
  end type yTruncatedCylinder

  type, public, extends(surface) :: zTruncatedCylinder
    private
    type(zPlane), dimension(:), pointer   :: zPlanes => null()
    type(zCylinder), pointer              :: cyl => null()  ! substituent cylinder
    real(defReal)                         :: a              ! the half-width separation of the planes
    real(defReal), dimension(3)           :: origin = [ZERO, ZERO, ZERO]
  contains
    procedure :: init => initZTruncCyl
    procedure :: evaluate => evaluateZTruncCyl
    procedure :: distanceToSurface => distanceToZTruncCyl
    procedure :: normalVector => normalZTruncCyl
    procedure :: reflectiveTransform => reflectiveTransformZTruncCyl
    procedure :: whichSurface => whichSurfaceZTruncCyl
    procedure :: setBoundaryConditions => setBoundaryConditionsZTruncCylinder
    procedure :: boundaryTransform => boundaryTransformZTruncCylinder
  end type zTruncatedCylinder

contains

!!
!! X-truncated cylinder procedures
!!
  !!
  !! Initialise the box as six plane surfaces
  !!
  subroutine initXTruncCyl(self, origin, a, radius, id, name)
    class(xTruncatedCylinder), intent(inout) :: self
    real(defReal), dimension(3), intent(in)  :: origin
    real(defReal), intent(in)                :: a, &
                                                radius
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

  end subroutine initXTruncCyl

  !!
  !! Evaluate the surface function of the square cylinder
  !!
  function evaluateXTruncCyl(self, r) result(res)
    class(xTruncatedCylinder), intent(in)   :: self
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

  end function evaluateXTruncCyl

  !!
  !! Calculate the distance to the nearest surface of the cylinder
  !! Requires checking that only real surfaces are intercepted,
  !! i.e., not the extensions of the plane or cylinderical surfaces
  !!
  function distanceToXTruncCyl(self, r, u)result(distance)
    class(xTruncatedCylinder), intent(in)   :: self
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

  end function distanceToXTruncCyl

  !!
  !! Apply a reflective transformation to a particle during delta tracking
  !! Do so by determining which plane the particle intersects and applying the plane reflection
  !!
  !! This routine is obviated due to the implementation in the transport operator and cell
  !!
  subroutine reflectiveTransformXTruncCyl(self, r, u)
    class(xTruncatedCylinder), intent(in)      :: self
    real(defReal), dimension(3), intent(inout) :: r, u
    class(surface), pointer                    :: surfPointer

    call fatalError('reflectiveTransformXTruncCyl','This routine should not be called')
    surfPointer => self % whichSurface(r, u)
    call surfPointer % reflectiveTransform(r,u)

  end subroutine reflectiveTransformXTruncCyl

  !!
  !! Determine on which surface the particle is located and obtain
  !! its normal vector
  !!
  function normalXTruncCyl(self, r) result(normal)
    class(xTruncatedCylinder), intent(in)   :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3)             :: normal
    real(defReal)                           :: posBound, negBound

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
      call fatalError('normalVector, xTruncatedCylinder','Point is not on a surface')
    end if

  end function normalXTruncCyl

  !!
  !! Helper routine: find the plane which a particle will have crossed
  !! This is done by calculating the distance to each surface
  !! Include inside or outside to allow generality (need to change
  !! particle direction otherwise)
  !!
  function whichSurfaceXTruncCyl(self, r, u) result(surfPointer)
    class(xTruncatedCylinder), intent(in)   :: self
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

  end function whichSurfaceXTruncCyl

  !!
  !! Set boundary conditions for an xTruncCylinder
  !!
  subroutine setBoundaryConditionsXTruncCylinder(self, BC)
    class(xTruncatedCylinder), intent(inout)    :: self
    integer(shortInt), dimension(6), intent(in) :: BC

    ! Positive x boundary
    if(BC(1) == vacuum) then
      self % xPlanes(1) % isVacuum = .TRUE.
    else if(BC(1) == reflective) then
      self % xPlanes(1) % isReflective = .TRUE.
    else if(BC(1) == periodic) then
      if(BC(2) /= periodic) then
        call fatalError('setBoundaryConditionsXTruncCylinder', &
        'Both positive and negative boundary conditions must be periodic')
      else
        self % xPlanes(1) % isPeriodic = .TRUE.
        self % xPlanes(1) % periodicTranslation = [-TWO*self % a, ZERO, ZERO]
      end if
    else
      call fatalError('setBoundaryConditionsXTruncCylinder','Invalid boundary condition provided')
    end if

    ! Negative x boundary
    if(BC(2) == vacuum) then
      self % xPlanes(2) % isVacuum = .TRUE.
    else if(BC(2) == reflective) then
      self % xPlanes(2) % isReflective = .TRUE.
    else if(BC(2) == periodic) then
      if(BC(1) /= periodic) then
        call fatalError('setBoundaryConditionsXTruncCylinder', &
        'Both positive and negative boundary conditions must be periodic')
      else
        self % xPlanes(2) % isPeriodic = .TRUE.
        self % xPlanes(2) % periodicTranslation = [TWO*self % a, ZERO, ZERO]
      end if
    else
      call fatalError('setBoundaryConditionsXTruncCylinder','Invalid boundary condition provided')
    end if

    if(any(BC(3:6) /= vacuum))then
      call fatalError('setBoundaryConditionsXTruncCylinder','Cylinder boundaries may only be vacuum')
    else
      self % cyl % isVacuum = .TRUE.
    end if
  end subroutine setBoundaryConditionsXTruncCylinder

  !!
  !! Apply boundary transformation
  !!
  subroutine boundaryTransformXTruncCylinder(self, r, u, isVacuum)
    class(xTruncatedCylinder), intent(in)      :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u
    logical(defBool), intent(inout)            :: isVacuum
    logical(defBool)                           :: front, back, outsideCyl

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
      call fatalError('boundaryTransform, xTruncatedCylinder',&
      'Cannot apply boundary condition: point is inside the surface')
    end if

  end subroutine boundaryTransformXTruncCylinder

!!
!! Y-truncated cylinder procedures
!!
  !!
  !! Initialise the box as six plane surfaces
  !!
  subroutine initYTruncCyl(self, origin, a, radius, id, name)
    class(yTruncatedCylinder), intent(inout) :: self
    real(defReal), dimension(3), intent(in)  :: origin
    real(defReal), intent(in)                :: a, &
                                                radius
    integer(shortInt), intent(in), optional  :: id
    character(*), optional, intent(in)       :: name

    self % isCompound = .true.
    self % a = a
    if(a < surface_tol) &
    call fatalError('init, yTruncCylinder','Width must be greater than surface tolerance')
    if(radius < surface_tol) &
    call fatalError('init, yTruncCylinder','Radius must be greater than surface tolerance')
    self % origin = origin
    if(present(id)) self % id = id
    if(present(name)) self % name = name

    if(associated(self % yPlanes)) deallocate (self % yPlanes)
    allocate(self%yPlanes(2))

    call self % yPlanes(1) % init(origin(2) + a)
    call self % yPlanes(2) % init(origin(2) - a)
    call self % cyl % init(radius, origin)

  end subroutine initYTruncCyl

  !!
  !! Evaluate the surface function of the square cylinder
  !!
  function evaluateYTruncCyl(self, r) result(res)
    class(yTruncatedCylinder), intent(in)   :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: res
    real(defReal), dimension(2)             :: resY      ! Results from yPlanes
    real(defReal)                           :: resCyl    ! Results from cylinder
    real(defReal)                           :: absMinRes ! Identify if particle sits on a surface
    real(defReal)                           :: testRes   ! The most positive residual

    ! Evaluate the top plane's surface function
    resY(1) = self % yPlanes(1) % evaluate(r)

    ! If the result is greater than the surface tolerance
    ! then the particle is not inside the cylinder and the result should
    ! be returned
    if (resY(1) > surface_tol) then
      res = resY(1)
      return
    end if

    resY(2) = -self % yPlanes(2) % evaluate(r)
    absMinRes = min(abs(resY(1)),abs(resY(2)))
    testRes = max(resY(1),resY(2))

    ! Same as before - particle is outside the cylinder
    if (resY(2) > surface_tol) then
      res = resY(2)
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

  end function evaluateYTruncCyl

  !!
  !! Calculate the distance to the nearest surface of the cylinder
  !! Requires checking that only real surfaces are intercepted,
  !! i.e., not the extensions of the plane or cylinderical surfaces
  !!
  function distanceToYTruncCyl(self, r, u)result(distance)
    class(yTruncatedCylinder), intent(in) :: self
    real(defReal), dimension(3), intent(in) :: r, u
    real(defReal), dimension(3) :: testPoint
    real(defReal) :: posBound, negBound
    real(defReal) :: distance
    real(defReal) :: testDistance
    real(defReal) :: testCyl

    distance = INFINITY
    ! Find the positive and negative bounds which the particle
    ! must fall within
    posBound = self % yPlanes(1) % y0
    negBound = self % yPlanes(2) % y0

    testDistance = self%yPlanes(1)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    testCyl = self % cyl % evaluate(testPoint)

    ! Ensure point is within the radius of the cylinder
    if (testCyl < surface_tol)  then
      if (testDistance < distance) distance = testDistance
    end if

    testDistance = self%yPlanes(2)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    testCyl = self % cyl % evaluate(testPoint)

    if (testCyl < surface_tol) then
      if (testDistance < distance) distance = testDistance
    end if

    testDistance = self% cyl %distanceToSurface(r,u)
    testPoint = r + u*testDistance
    if ((testPoint(2) < posBound) .and. (testPoint(2) > negBound)) then
      if (testDistance < distance) distance = testDistance
    end if

  end function distanceToYTruncCyl

  !
  !! Apply a reflective transformation to a particle during delta tracking
  !! Do so by determining which plane the particle intersects and applying the plane reflection
  !!
  !! This routine is obviated due to the implementation in the transport operator and cell
  !!
  subroutine reflectiveTransformYTruncCyl(self, r, u)
    class(yTruncatedCylinder), intent(in)      :: self
    real(defReal), dimension(3), intent(inout) :: r, u
    class(surface), pointer                    :: surfPointer

    call fatalError('reflectiveTransformYTruncCyl','This routine should not be called')
    surfPointer => self % whichSurface(r, u)
    call surfPointer % reflectiveTransform(r,u)

  end subroutine reflectiveTransformYTruncCyl

  !!
  !! Determine on which surface the particle is located and obtain
  !! its normal vector
  !!
  function normalYTruncCyl(self, r)result(normal)
    class(yTruncatedCylinder), intent(in)   :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3)             :: normal
    real(defReal)                           :: posBound, negBound

    ! Compare the point's position to the maximum and minimum
    posBound = self % yPlanes(1) % y0
    negBound = self % yPlanes(2) % y0

    if (abs(posBound - r(2)) < surface_tol) then
      normal = self % yPlanes(1) % normalVector(r)
      return
    else if (abs(negBound - r(2)) < surface_tol) then
      normal = self % yPlanes(2) % normalVector(r)
      return
    else if (abs(self % cyl % evaluate(r)) < surface_tol) then
      normal = self % cyl % normalVector(r)
      return
    else
      call fatalError('normalVector, yTruncatedCylinder','Point is not on a surface')
    end if

  end function normalYTruncCyl

  !!
  !! Helper routine: find the plane which a particle will have crossed
  !! This is done by calculating the distance to each surface
  !! Include inside or outside to allow generality (need to change
  !! particle direction otherwise)
  !!
  function whichSurfaceYTruncCyl(self, r, u) result(surfPointer)
    class(yTruncatedCylinder), intent(in)   :: self
    real(defReal), dimension(3), intent(in) :: r, u
    class(surface), pointer                 :: surfPointer
    real(defReal), dimension(3)             :: testPoint
    real(defReal)                           :: posBound, negBound
    real(defReal)                           :: testDistance, distance
    real(defReal)                           :: testCyl

    distance = INFINITY
    posBound = self % yPlanes(1) % y0
    negBound = self % yPlanes(2) % y0

    ! Evaluate distance to each plane and point to surface
    ! with the minimum real distance
    testDistance = self%yPlanes(1)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    testCyl = self % cyl % evaluate(testPoint)

    if (testCyl < surface_tol) then
      if (testDistance < distance) then
        distance = testDistance
        surfPointer => self % yPlanes(1)
      end if
    end if

    testDistance = self%yPlanes(2)%distanceToSurface(r,u)
    testPoint = r + u*testDistance
    testCyl = self % cyl % evaluate(testPoint)

    if (testCyl < surface_tol) then
      if (testDistance < distance) then
        distance = testDistance
        surfPointer => self % yPlanes(2)
      end if
    end if

    testDistance = self % cyl % distanceToSurface(r,u)
    testPoint = r + u*testDistance
    if ((testPoint(2) < posBound).and.(testPoint(2) > negBound))  then
      if (testDistance < distance) then
        distance = testDistance
        surfPointer => self % cyl
      end if
    end if

  end function whichSurfaceYTruncCyl

  !!
  !! Set boundary conditions for a yTruncCylinder
  !!
  subroutine setBoundaryConditionsYTruncCylinder(self, BC)
    class(yTruncatedCylinder), intent(inout)    :: self
    integer(shortInt), dimension(6), intent(in) :: BC

    ! Positive y boundary
    if(BC(3) == vacuum) then
      self % yPlanes(1) % isVacuum = .TRUE.
    else if(BC(3) == reflective) then
      self % yPlanes(1) % isReflective = .TRUE.
    else if(BC(3) == periodic) then
      if(BC(4) /= periodic) then
        call fatalError('setBoundaryConditionsYTruncCylinder', &
        'Both positive and negative boundary conditions must be periodic')
      else
        self % yPlanes(1) % isPeriodic = .TRUE.
        self % yPlanes(1) % periodicTranslation = [ZERO, -TWO*self % a, ZERO]
      end if
    else
      call fatalError('setBoundaryConditionsYTruncCylinder','Invalid boundary condition provided')
    end if

    ! Negative y boundary
    if(BC(4) == vacuum) then
      self % yPlanes(2) % isVacuum = .TRUE.
    else if(BC(4) == reflective) then
      self % yPlanes(2) % isReflective = .TRUE.
    else if(BC(4) == periodic) then
      if(BC(3) /= periodic) then
        call fatalError('setBoundaryConditionsYTruncCylinder', &
        'Both positive and negative boundary conditions must be periodic')
      else
        self % yPlanes(2) % isPeriodic = .TRUE.
        self % yPlanes(2) % periodicTranslation = [ZERO, TWO*self % a, ZERO]
      end if
    else
      call fatalError('setBoundaryConditionsYTruncCylinder','Invalid boundary condition provided')
    end if

    if(any(BC(1:2) /= vacuum) .OR. any(BC(5:6) /= vacuum))then
      call fatalError('setBoundaryConditionsYTruncCylinder','Cylinder boundaries may only be vacuum')
    else
      self % cyl % isVacuum = .TRUE.
    end if
  end subroutine setBoundaryConditionsYTruncCylinder

  !!
  !! Apply boundary transformation
  !!
  subroutine boundaryTransformYTruncCylinder(self, r, u, isVacuum)
    class(yTruncatedCylinder), intent(in)      :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u
    logical(defBool), intent(inout)            :: isVacuum
    logical(defBool)                           :: left, right, outsideCyl

    outsideCyl = self % cyl % halfspace(r, u)
    if (outsideCyl) then
      call self % cyl % boundaryTransform(r, u, isVacuum)
      return
    end if

    left = self % yPlanes(1) % halfspace(r, u)
    right = .NOT. self % yPlanes(2) % halfspace(r, u)
    if (left) then
      call self % yPlanes(1) % boundaryTransform(r, u, isVacuum)
    else if (right) then
      call self % yPlanes(2) % boundaryTransform(r, u, isVacuum)
    else
      call fatalError('boundaryTransform, yTruncatedCylinder',&
      'Cannot apply boundary condition: point is inside the surface')
    end if

  end subroutine boundaryTransformYTruncCylinder

!!
!! Z-truncated cylinder procedures
!!
  !!
  !! Initialise the box as six plane surfaces
  !!
  subroutine initZTruncCyl(self, origin, a, radius, id, name)
    class(zTruncatedCylinder), intent(inout) :: self
    real(defReal), dimension(3), intent(in)  :: origin
    real(defReal), intent(in)                :: a, &
                                                radius
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

  end subroutine initZTruncCyl

  !!
  !! Evaluate the surface function of the square cylinder
  !!
  function evaluateZTruncCyl(self, r) result(res)
    class(zTruncatedCylinder), intent(in)   :: self
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

  end function evaluateZTruncCyl

  !!
  !! Calculate the distance to the nearest surface of the cylinder
  !! Requires checking that only real surfaces are intercepted,
  !! i.e., not the extensions of the plane or cylinderical surfaces
  !!
  function distanceToZTruncCyl(self, r, u) result(distance)
    class(zTruncatedCylinder), intent(in)   :: self
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

  end function distanceToZTruncCyl

  !!
  !! Apply a reflective transformation to a particle during delta tracking
  !! Do so by determining which plane the particle intersects and applying the plane reflection
  !!
  !! This route is obviated due to the implementation in transport operator and cell
  !!
  subroutine reflectiveTransformZTruncCyl(self, r, u)
    class(zTruncatedCylinder), intent(in)      :: self
    real(defReal), dimension(3), intent(inout) :: r, u
    class(surface), pointer                    :: surfPointer

    call fatalError('reflectiveTransformZTruncCyl','This routine should not be called')
    surfPointer => self % whichSurface(r, u)
    call surfPointer % reflectiveTransform(r,u)

  end subroutine reflectiveTransformZTruncCyl

  !!
  !! Determine on which surface the particle is located and obtain
  !! its normal vector
  !!
  function normalZTruncCyl(self, r) result(normal)
    class(zTruncatedCylinder), intent(in)   :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3)             :: normal
    real(defReal)                           :: posBound, negBound

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
      call fatalError('normalVector, zTruncatedCylinder','Point is not on a surface')
    end if

  end function normalZTruncCyl

  !!
  !! Helper routine: find the plane which a particle will have crossed
  !! This is done by calculating the distance to each surface
  !! Include inside or outside to allow generality (need to change
  !! particle direction otherwise)
  !!
  function whichSurfaceZTruncCyl(self, r, u) result(surfPointer)
    class(zTruncatedCylinder), intent(in)   :: self
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

  end function whichSurfaceZTruncCyl

  !!
  !! Set boundary conditions for a zTruncCylinder
  !!
  subroutine setBoundaryConditionsZTruncCylinder(self, BC)
    class(zTruncatedCylinder), intent(inout) :: self
    integer(shortInt), dimension(6), intent(in) :: BC

    ! Positive z boundary
    if(BC(5) == vacuum) then
      self % zPlanes(1) % isVacuum = .TRUE.
    else if(BC(5) == reflective) then
      self % zPlanes(1) % isReflective = .TRUE.
    else if(BC(5) == periodic) then
      if(BC(6) /= periodic) then
        call fatalError('setBoundaryConditionsZTruncCylinder', &
        'Both positive and negative boundary conditions must be periodic')
      else
        self % zPlanes(1) % isPeriodic = .TRUE.
        self % zPlanes(1) % periodicTranslation = [ZERO, ZERO, -TWO*self % a]
      end if
    else
      call fatalError('setBoundaryConditionsZTruncCylinder','Invalid boundary condition provided')
    end if

    ! Negative z boundary
    if(BC(6) == vacuum) then
      self % zPlanes(2) % isVacuum = .TRUE.
    else if(BC(6) == reflective) then
      self % zPlanes(2) % isReflective = .TRUE.
    else if(BC(6) == periodic) then
      if(BC(5) /= periodic) then
        call fatalError('setBoundaryConditionsZTruncCylinder', &
        'Both positive and negative boundary conditions must be periodic')
      else
        self % zPlanes(2) % isPeriodic = .TRUE.
        self % zPlanes(2) % periodicTranslation = [ZERO, ZERO, TWO*self % a]
      end if
    else
      call fatalError('setBoundaryConditionsZTruncCylinder','Invalid boundary condition provided')
    end if

    if(any(BC(1:4) /= vacuum))then
      call fatalError('setBoundaryConditionsZTruncCylinder','Cylinder boundaries may only be vacuum')
    else
      self % cyl % isVacuum = .TRUE.
    end if
  end subroutine setBoundaryConditionsZTruncCylinder

  !!
  !! Apply boundary transformation
  !!
  subroutine boundaryTransformZTruncCylinder(self, r, u, isVacuum)
    class(zTruncatedCylinder), intent(in)      :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u
    logical(defBool), intent(inout)            :: isVacuum
    logical(defBool)                           :: above, below, outsideCyl

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
      call fatalError('boundaryTransform, zTruncatedCylinder',&
      'Cannot apply boundary condition: point is inside the surface')
    end if

  end subroutine boundaryTransformZTruncCylinder

end module truncatedCylinder_class
