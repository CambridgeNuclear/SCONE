module infSurf_class

  use numPrecision
  use universalVariables
  use genericProcedures, only : fatalError
  use dictionary_class,  only : dictionary
  use surface_inter,     only : surface, printSurfDef

  implicit none
  private

  !!
  !! Constants describing surface properties
  !!
  character(nameLen),parameter :: TYPE_NAME    = 'infSurf'

  !!
  !! Constructor
  !!
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
    procedure :: type
    procedure :: getDef
    procedure :: distanceToSurface
    procedure :: normalVector
    procedure :: whichSurface
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
  !! Return parameter character containing TYPE NAME
  !!
  function type(self)
    class(infSurf), intent(in) :: self
    character(nameLen)         :: type

    type = TYPE_NAME

  end function type

  !!
  !! Return string with definition of this surface
  !!
  pure subroutine getDef(self,string)
    class(infSurf), intent(in)             :: self
    character(:),allocatable,intent(inout) :: string

    string = printSurfDef(TYPE_NAME, [ZERO])

  end subroutine getDef

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
