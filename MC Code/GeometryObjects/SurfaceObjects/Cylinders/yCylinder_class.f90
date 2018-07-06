module yCylinder_class

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
  character(nameLen),parameter :: TYPE_NAME    = 'yCylinder'
  logical(defBool),parameter   :: ACCEPTS_BC   = .true.

  !!
  !! Constructor
  !!
  interface yCylinder
    module procedure yCylinder_fromDict
  end interface

  !!
  !! Cylinder aligned with Y-axis
  !!
  type, public, extends(surface) :: yCylinder
    real(defReal)                :: rSquared = ZERO             ! squared radius of the cylinder
    real(defReal)                :: radius = ZERO               ! radius of the cylinder
    real(defReal), dimension(3)  :: origin = [ZERO, ZERO, ZERO] ! 2D origin of cylinder assumed parallel to y-axis

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

  end type yCylinder

contains
  !!
  !! Initialise Y-Cylinder from components
  !!
  subroutine init(self, radius, origin, id, name)
    class(yCylinder), intent(inout)         :: self
    real(defReal), intent(in)               :: radius
    real(defReal), dimension(3), intent(in) :: origin
    integer(shortInt), intent(in), optional :: id
    character(*), optional, intent(in)      :: name

    self % radius = radius
    if (radius < surface_tol) &
    call fatalError('init, yCylinder','Radius must be greater than surface tolerance')
    self % rSquared = radius*radius
    self % origin = origin
    if(present(id)) self % id = id
    if(present(name)) self % name = name

    ! Hash and store surface definition
    call self % hashSurfDef()

  end subroutine init

  !!
  !! Returns and initialised instance of yCylinder from dictionary and name
  !!
  function yCylinder_fromDict(dict,name) result(new)
    class(dictionary), intent(in)  :: dict
    character(nameLen), intent(in) :: name
    type(yCylinder)                :: new
    integer(shortInt)              :: id
    real(defReal)                  :: radius
    real(defReal), dimension(3)    :: origin
    character(100),parameter :: Here ='yCylinder_fromDict ( yCylinder_class.f90)'

    id = dict % getInt('id')
    if(id < 1) call fatalError(Here,'Invalid surface id provided')

    radius = dict % getReal('radius')
    origin = dict % getRealArray('origin')
    call new % init(radius, origin, id, name)

  end function yCylinder_fromDict

  !!
  !! Evaluate remainder of cylinder equation
  !!
  function evaluate(self, r) result(res)
    class(yCylinder), intent(in)            :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: res

    res = (r(1) - self%origin(1))**2 + (r(3) - self%origin(3))**2 - self%rSquared

  end function evaluate

  !!
  !! Return parameter character containing TYPE NAME
  !!
  function type(self)
    class(yCylinder), intent(in) :: self
    character(nameLen)           :: type

    type = TYPE_NAME

  end function type

  !!
  !! Return string with definition of this surface
  !!
  pure subroutine getDef(self,string)
    class(yCylinder), intent(in)           :: self
    character(:),allocatable,intent(inout) :: string

    string = printSurfDef(TYPE_NAME, [self % radius, self % origin])

  end subroutine getDef

  !!
  !! Override base type function to returns .false.
  !!
  function cannotBeBoundary(self) result(itCant)
    class(yCylinder), intent(in) :: self
    logical(defBool)             :: itCant

    itCant = .not.ACCEPTS_BC

  end function cannotBeBoundary

  !!
  !! Calculate distance to cylinder along direction u
  !!
  function distanceToSurface(self, r, u) result(distance)
    class(yCylinder), intent(in)            :: self
    real(defReal), dimension(3), intent(in) :: r, u
    real(defReal)                           :: xBar, zBar
    real(defReal)                           :: k, c, a
    real(defReal)                           :: discriminant, distance

    xBar = r(1) - self%origin(1)
    zBar = r(3) - self%origin(3)

    k = xBar*u(1) + zBar*u(3)
    a = ONE - u(2)*u(2) ! = u(1)*u(1) + u(3)*u(3)
    c = xBar*xBar + zBar*zBar - self%rSquared
    discriminant = k*k - a*c

    ! No intersection
    if ((a==ZERO).OR.(discriminant<ZERO)) then
      distance = INFINITY
      return
    ! Particle is on the cylinder - the distance will be either positive or negative
    else if (abs(c) < surface_tol) then
      if (k >= ZERO) then
        distance = INFINITY
        return
      else
        distance = (-k + sqrt(discriminant))/a
        return
      end if
    ! Particle is inside - take solution with + before sqrt
    else if (c<ZERO) then
      distance = (-k + sqrt(discriminant))/a
      return
    ! Particle is outside - both distances are either positive or negative
    else
      distance = (-k -sqrt(discriminant))/a
      if (distance<ZERO) distance = INFINITY
      return
    end if

  end function distanceToSurface

  !!
  !! Return normal to the cylinder
  !!
  function normalVector(self,r) result(normal)
    class(yCylinder), intent(in)            :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3)             :: normal

    normal = [TWO*(r(1) - self%origin(1)), ZERO, TWO*(r(3) - self%origin(3))]

  end function normalVector

  !!
  !! Give an error: this routine should not be called for a non-compound surface
  !!
  function whichSurface(self, r, u) result(surfPointer)
    class(yCylinder), intent(in)            :: self
    real(defReal), dimension(3), intent(in) :: r, u
    class(surface), pointer                 :: surfPointer
    character(100),parameter :: Here ='whichSurface (yCylinder_class.f90)'

    call fatalError(Here,'This function should never be called for a simple surface')

  end function whichSurface

  !!
  !! Set boundary conditions for an yCylinder: may only be vacuum
  !!
  subroutine setBoundaryConditions(self, BC)
    class(yCylinder), intent(inout)             :: self
    integer(shortInt), dimension(:), intent(in) :: BC
    character(100),parameter :: Here ='setBoundaryConditions (yCylinder_class.f90)'

    if (any(BC /= vacuum)) then
      call fatalError(Here,'Cylinder boundaries may only be vacuum')
    else
      self % isVacuum = .TRUE.
    end if

  end subroutine setBoundaryConditions

  !!
  !! Apply boundary transformation
  !!
  subroutine boundaryTransform(self, r, u, isVacuum)
    class(yCylinder), intent(in) :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u
    logical(defBool), intent(inout)            :: isVacuum
    character(100),parameter :: Here ='boundaryTransform (yCylinder_class.f90)'

    if (self % isVacuum) then
      isVacuum = .TRUE.

    else
      call fatalError(Here,'This should only be called for a cylinder with vacuum boundaries')

    end if

  end subroutine boundaryTransform

    
end module yCylinder_class
