module zPlane_class

  use numPrecision
  use universalVariables
  use genericProcedures,    only : fatalError, dotProduct

  use surface_inter,        only : surface

  implicit none
  private

  !!
  !! Simple plane perpendicular to z-axis
  !!
  type, public, extends(surface) :: zPlane
    real(defReal) :: z0 ! z-axis offset

  contains
    procedure :: init
    procedure :: evaluate
    procedure :: distanceToSurface
    procedure :: reflectiveTransform
    procedure :: normalVector
    procedure :: whichSurface
    procedure :: setBoundaryConditions
    procedure :: boundaryTransform

  end type zPlane

contains

  !!
  !! Initialise Z-Plane from offest
  !!
  subroutine init(self, z0, id, name)
    class(zPlane), intent(inout)            :: self
    real(defReal), intent(in)               :: z0
    integer(shortInt), intent(in), optional :: id
    character(*), optional, intent(in)      :: name

    self % z0 = z0
    if(present(id)) self % id = id
    if(present(name)) self % name = name

  end subroutine init

  !!
  !! Evaluate plane distance
  !!
  function evaluate(self, r) result(res)
    class(zPlane), intent(in)               :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: res

    res = r(3) - self % z0

  end function evaluate

  !!
  !! Calculate distance to plane along direction u
  !!
  function distanceToSurface(self,r,u)result(distance)
    class(zPlane), intent(in)               :: self
    real(defReal), dimension(3), intent(in) :: r, u
    real(defReal)                           :: distance, w

    w = u(3)
    distance = self%z0 - r(3)
    if ((w==ZERO) .OR. (abs(distance) < surface_tol)) then
      distance = INFINITY
      return
    end if

    distance = (self % z0 - r(3))/w
    if (distance < ZERO) distance = INFINITY
    return

  end function distanceToSurface

  !!
  !! Perform reflection
  !!
  subroutine reflectiveTransform(self,r,u)
    class(zPlane), intent(in)                  :: self
    real(defReal), dimension(3), intent(inout) :: r, u
    real(defReal) :: zDisplacement

    ! Reflect the particle z-coordinate across the plane
    zDisplacement = r(3) - self%z0
    r(3) = r(3) - TWO*zDisplacement

    ! Reflect the particle direction (independent of intersection point for plane)
    u(3) = -u(3)

  end subroutine reflectiveTransform

  !!
  !! Returns vector normal to the plane
  !!
  function normalVector(self,r)result(normal)
    class(zPlane), intent(in)               :: self
    real(defReal), dimension(3)             :: normal
    real(defReal), dimension(3), intent(in) :: r

    normal = [ZERO, ZERO, ONE]

  end function normalVector

  !!
  !! Give an error: this routine should not be called for a non-compound surface
  !!
  function whichSurface(self, r, u) result(surfPointer)
    class(zPlane), intent(in)               :: self
    real(defReal), dimension(3), intent(in) :: r, u
    class(surface), pointer                 :: surfPointer

    call fatalError('whichZPlane','This function should never be called for a simple surface')

  end function whichSurface

  !!
  !! Crash on attempting to set boundary conditions for a plane object
  !!
  subroutine setBoundaryConditions(self, BC)
    class(zPlane), intent(inout)                :: self
    integer(shortInt), dimension(6), intent(in) :: BC

    call fatalError('setBoundaryConditionsZPlane','Boundary conditions may not be set for a plane surface')

  end subroutine setBoundaryConditions

  !!
  !! Apply boundary transformations
  !!
  subroutine boundaryTransform(self, r, u, isVacuum)
    class(zPlane), intent(in)                  :: self
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
      call fatalError('boundaryTransform, zPlane','No boundary condition applied to surface')
    end if

  end subroutine boundaryTransform

end module zPlane_class
