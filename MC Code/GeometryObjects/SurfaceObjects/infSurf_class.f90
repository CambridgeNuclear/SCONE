module infSurf_class

  use numPrecision
  use universalVariables
  use genericProcedures, only : fatalError
  use dictionary_class,  only : dictionary
  use surface_inter,     only : surface

  implicit none
  private

  interface infSurf
    module procedure infSurf_fromDict
  end interface

  !!
  !! Infinite surface for homogeneous regions
  !!
  type, public, extends (surface) :: infSurf
    private
  contains
    procedure :: init
    procedure :: evaluate
    procedure :: distanceToSurface
    procedure :: reflectiveTransform
    procedure :: normalVector
    procedure :: whichSurface
    procedure :: setBoundaryConditions
    procedure :: boundaryTransform
  end type infSurf

contains

  !!
  !! Create an infinity surface object
  !!
  subroutine init(self, id, name)
    class(infSurf), intent(inout)           :: self
    integer(shortInt), intent(in), optional :: id
    character(*), optional, intent(in)      :: name

    if(present(id))   self % id = id
    if(present(name)) self % name = name

  end subroutine init

  !!
  !! Returns an initialised instance of infSurf from dictionary and name
  !!
  function infSurf_fromDict(dict,name) result(new)
    class(dictionary), intent(in)  :: dict
    character(nameLen),intent(in)  :: name
    type(infSurf)                  :: new
    integer(shortInt)              :: id
    character(100), parameter :: Here ='infSurf_fromDict (infSurf_class.f90)'

    id = dict % getInt('id')
    if(id < 1) call fatalError(Here,'Invalid surface id provided')

    call new % init(id,name)

  end function infSurf_fromDict

  !!
  !! Always return inside for an infinite surface
  !!
  function evaluate(self, r) result(res)
    class(infSurf), intent(in)              :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: res

    res = -ONE

  end function evaluate

  !!
  !! Calculate distance to infinity's surface - always infinity
  !!
  function distanceToSurface(self,r,u) result(distance)
    class(infSurf), intent(in)              :: self
    real(defReal), dimension(3), intent(in) :: r, u
    real(defReal)                           :: distance

    distance = INFINITY

  end function distanceToSurface

  !!
  !! Perform a co-ordinate transform on a particle to apply reflective boundary condition
  !!
  subroutine reflectiveTransform(self, r, u)
    class(infSurf), intent(in)                 :: self
    real(defReal), dimension(3), intent(inout) :: r, u
    character(100), parameter :: Here ='reflectiveTransform (infSurf_class.f90)'

    ! Reflective transforms will not be allowed to occur in geometries other than planes
    call fatalError(Here,'Infinite surfaces may not have reflective boundaries')

  end subroutine reflectiveTransform

  !!
  !! Supply the normal vector
  !!
  function normalVector(self, r) result(normal)
    class(infSurf), intent(in)              :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3)             :: normal
    character(100), parameter :: Here ='normalVector (infSurf_class.f90)'

    call fatalError(Here,'Infinite surfaces do not have normals')

  end function normalVector

  !!
  !! Give an error: this routine should not be called for a non-compound surface
  !!
  function whichSurface(self, r, u) result(surfPointer)
    class(infSurf), intent(in)              :: self
    real(defReal), dimension(3), intent(in) :: r, u
    class(surface), pointer                 :: surfPointer
    character(100), parameter :: Here ='whichSurface (infSurf_class.f90)'

    call fatalError(Here,'This function should never be called for a simple surface')

  end function whichSurface

  !!
  !! Set boundary conditions for an infinite surface: may only be vacuum
  !! Doesn't matter - a particle will never reach the surface anyway
  !!
  subroutine setBoundaryConditions(self, BC)
    class(infSurf), intent(inout)               :: self
    integer(shortInt), dimension(6), intent(in) :: BC
    character(100), parameter :: Here ='setBoundaryConditions (infSurf_class.f90)'

    if (any(BC /= vacuum)) then
      call fatalError(Here,'Infinite surfaces may only be vacuum')

    else
      self % isVacuum = .TRUE.

    end if
  end subroutine setBoundaryConditions

  !!
  !! Apply boundary transformation
  !!
  subroutine boundaryTransform(self, r, u, isVacuum)
    class(infSurf), intent(in)                 :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u
    logical(defBool), intent(inout)            :: isVacuum
    character(100), parameter :: Here ='boundaryTransform (infSurf_class.f90)'

    call fatalError(Here,'Infinite surfaces should not have associated boundary conditions')

  end subroutine boundaryTransform

end module infSurf_class
