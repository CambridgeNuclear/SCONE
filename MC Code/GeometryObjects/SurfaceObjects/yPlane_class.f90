module yPlane_class

  use numPrecision
  use universalVariables
  use genericProcedures,    only : fatalError, dotProduct

  use surface_inter,        only : surface

  implicit none
  private

  !!
  !! Simple plane perpendicular to y-axis
  !!
  type, public, extends(surface) :: yPlane
    real(defReal) :: y0 ! y-axis offset

  contains
    procedure :: init
    procedure :: evaluate
    procedure :: distanceToSurface
    procedure :: reflectiveTransform
    procedure :: normalVector
    procedure :: whichSurface
    procedure :: setBoundaryConditions
    procedure :: boundaryTransform

  end type yPlane


contains

  !!
  !! Initialise Y-Plabne from offset
  !!
  subroutine init(self, y0, id, name)
    class(yPlane), intent(inout)            :: self
    real(defReal), intent(in)               :: y0
    integer(shortInt), intent(in), optional :: id
    character(*), optional, intent(in)      :: name

    self % y0 = y0
    if(present(id)) self % id = id
    if(present(name)) self % name = name

  end subroutine init

  !!
  !! Evaluate plane distance
  !!
  function evaluate(self, r) result(res)
    class(yPlane), intent(in)               :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: res

    res = r(2) - self % y0

  end function evaluate

  !!
  !! Calculate distance to plane along direction u
  !!
  function distanceToSurface(self,r,u)result(distance)
    class(yPlane), intent(in)               :: self
    real(defReal), dimension(3), intent(in) :: r, u
    real(defReal)                           :: distance, v

    v = u(2)
    distance = self%y0 - r(2)
    if ((v==ZERO) .OR. (abs(distance) < surface_tol)) then
      distance = INFINITY
      return
    end if

    distance = distance/v
    if (distance < ZERO) distance = INFINITY
    return

  end function distanceToSurface

  !!
  !! Perform reflection
  !!
  subroutine reflectiveTransform(self,r,u)
    class(yPlane), intent(in)                  :: self
    real(defReal), dimension(3), intent(inout) :: r, &
                                                  u
    real(defReal) :: yDisplacement

    ! Reflect the particle y-coordinate across the plane
    yDisplacement = r(2) - self%y0
    r(2) = r(2) - TWO*yDisplacement

    ! Reflect the particle direction (independent of intersection point for plane)
    u(2) = -u(2)

  end subroutine reflectiveTransform

  !!
  !! Returns vector normal to the plane
  !!
  function normalVector(self,r)result(normal)
    class(yPlane), intent(in)               :: self
    real(defReal), dimension(3)             :: normal
    real(defReal), dimension(3), intent(in) :: r

    normal = [ZERO, ONE, ZERO]
    return

  end function normalVector
  !!
  !! Give an error: this routine should not be called for a non-compound surface
  !!
  function whichSurface(self, r, u) result(surfPointer)
    class(yPlane), intent(in)               :: self
    real(defReal), dimension(3), intent(in) :: r, u
    class(surface), pointer                 :: surfPointer

    call fatalError('whichYPlane','This function should never be called for a simple surface')

  end function whichSurface

  !!
  !! Crash on attempting to set boundary conditions for a plane object
  !!
  subroutine setBoundaryConditions(self, BC)
    class(yPlane), intent(inout)                :: self
    integer(shortInt), dimension(6), intent(in) :: BC
    call fatalError('setBoundaryConditionsYPlane','Boundary conditions may not be set for a plane surface')
  end subroutine setBoundaryConditions

  !!
  !! Apply boundary transformations
  !!
  subroutine boundaryTransform(self, r, u, isVacuum)
    class(yPlane), intent(in)                  :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u
    logical(defBool), intent(inout)            :: isVacuum

    if (self % isVacuum) then
      isVacuum = .TRUE.
    else if (self % isPeriodic) then
      r = r + self % periodicTranslation
    else if (self % isReflective) then
      call self % reflectiveTransform(r,u)
    else
      call fatalError('boundaryTransform, yPlane','No boundary condition applied to surface')
    end if

  end subroutine boundaryTransform


end module yPlane_class
