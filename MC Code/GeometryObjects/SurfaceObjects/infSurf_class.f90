module infSurf_class

  use numPrecision
  use genericProcedures, only : fatalError
  use universalVariables

  use surface_class

  implicit none
  private

  !
  ! Infinite surface for homogeneous regions
  !
  type, public, extends (surface) :: infSurf
    private
  contains
    procedure :: init => initInf
    procedure :: evaluate => evaluateInf
    procedure :: distanceToSurface => distanceToInf
    procedure :: reflectiveTransform => reflectiveTransformInf
    procedure :: normalVector => normalVectorInf
    procedure :: whichSurface
    procedure :: setBoundaryConditions => setBoundaryConditionsInf
  end type infSurf

contains

  !!
  !! Create an infinity surface object
  !!
  subroutine initInf(self, id, name)
    class(infSurf), intent(inout)           :: self
    integer(shortInt), intent(in), optional :: id
    character(*), optional, intent(in)      :: name
    if(present(id)) self % id = id
    if(present(name)) self % name = name
  end subroutine initInf

  !!
  !! Always return inside for an infinite surface
  !!
  function evaluateInf(self, r) result(res)
    class(infSurf), intent(in)              :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: res
    res = -ONE
  end function evaluateInf

  !!
  !! Calculate distance to infinity's surface - always infinity
  !!
  function distanceToInf(self,r,u) result(distance)
    class(infSurf), intent(in)              :: self
    real(defReal), dimension(3), intent(in) :: r, u
    real(defReal)                           :: distance
    distance = INFINITY
  end function distanceToInf

  !!
  !! Perform a co-ordinate transform on a particle to apply reflective boundary condition
  !!
  subroutine reflectiveTransformInf(self, r, u)
    class(infSurf), intent(in)                 :: self
    real(defReal), dimension(3), intent(inout) :: r, u
    ! Reflective transforms will not be allowed to occur in geometries other than planes
    call fatalError('reflectiveTransform, infSurf','Infinite surfaces may not have reflective boundaries')
  end subroutine reflectiveTransformInf

  !!
  !! Supply the normal vector
  !!
  function normalVectorInf(self, r) result(normal)
    class(infSurf), intent(in)              :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3)             :: normal
    call fatalError('normalVectorInfSurf','Infinite surfaces do not have normals')
  end function normalVectorInf

  !!
  !! Give an error: this routine should not be called for a non-compound surface
  !!
  function whichSurface(self, r, u) result(surfPointer)
    class(infSurf), intent(in)              :: self
    real(defReal), dimension(3), intent(in) :: r, u
    class(surface), pointer                 :: surfPointer
    call fatalError('whichSurface, infSurf','This function should never be called for a simple surface')
  end function whichSurface

  !!
  !! Set boundary conditions for an infinite surface: may only be vacuum
  !! Doesn't matter - a particle will never reach the surface anyway
  !!
  subroutine setBoundaryConditionsInf(self, BC)
    class(infSurf), intent(inout)               :: self
    integer(shortInt), dimension(6), intent(in) :: BC

    if (any(BC /= vacuum)) then
      call fatalError('setBoundaryConditionsInf','Infinite surfaces may only be vacuum')
    else
      self % isVacuum = .TRUE.
    end if
  end subroutine setBoundaryConditionsInf

end module infSurf_class
