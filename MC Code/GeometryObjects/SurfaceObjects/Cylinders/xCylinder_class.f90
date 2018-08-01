module xCylinder_class

  use numPrecision
  use universalVariables
  use genericProcedures, only : fatalError
  use vector_class,      only : vector
  use dictionary_class,  only : dictionary
  use surface_inter,     only : surface, printSurfDef, surfaceShelf

  implicit none
  private

  !!
  !! Constants describing surface properties
  !!
  character(nameLen),parameter :: TYPE_NAME    = 'xCylinder'
  logical(defBool),parameter   :: ACCEPTS_BC   = .true.

  !!
  !! Constructor
  !!
  interface xCylinder
    module procedure xCylinder_fromDict
  end interface

  !!
  !! Cylinder aligned with x-axis
  !!
  type, public, extends(surface) :: xCylinder
    real(defReal)                :: rSquared = ZERO             ! squared radius of the cylinder
    real(defReal)                :: radius = ZERO               ! radius of the cylinder
    real(defReal), dimension(3)  :: origin = [ZERO, ZERO, ZERO] ! 2D origin of cylinder assumed parallel to x-axis
    integer(shortInt)            :: BC    = noBC
  contains
    ! Initialisation & Indentification procedures
    procedure :: init
    procedure :: type
    procedure :: getDef
    procedure :: boundingBox
    procedure :: cannotBeBoundary
    procedure :: setBoundaryConditions

    ! Runtime procedures
    procedure :: evaluate
    procedure :: distance
    procedure :: normalVector
    procedure :: boundaryTransform

  end type xCylinder

contains

  !!
  !! Initialise X-Cylinder from components
  !!
  subroutine init(self, radius, origin, id)
    class(xCylinder), intent(inout)         :: self
    real(defReal), intent(in)               :: radius
    real(defReal), dimension(3), intent(in) :: origin
    integer(shortInt), intent(in)           :: id

    self % radius = radius
    if(radius < surface_tol) &
    call fatalError('init, xCylinder','Radius must be greater than surface tolerance')
    self % rSquared = radius*radius
    self % origin = origin

    call self % setId(id)

    ! Hash and store surface definition
    call self % hashSurfDef()

  end subroutine init

  !!
  !! Returns and initialised instance of xCylinder from dictionary and name
  !!
  function xCylinder_fromDict(dict) result(new)
    class(dictionary), intent(in)  :: dict
    type(xCylinder)                :: new
    integer(shortInt)              :: id
    real(defReal)                  :: radius
    real(defReal), dimension(3)    :: origin
    character(100),parameter :: Here ='xCylinder_fromDict ( xCylinder_class.f90)'

    id = dict % getInt('id')
    if(id < 1) call fatalError(Here,'Invalid surface id provided')

    radius = dict % getReal('radius')
    origin = dict % getRealArray('origin')
    call new % init(radius, origin, id)

  end function xCylinder_fromDict

  !!
  !! Evaluate remainder of cylinder equation
  !!
  elemental subroutine evaluate(self, res, r)
    class(xCylinder), intent(in)   :: self
    real(defReal), intent(out)     :: res
    type(vector), intent(in)       :: r

    res = (r % v(2) - self % origin(2))**2 + (r % v(3) - self%origin(3))**2 - self % rSquared

  end subroutine evaluate

  !!
  !! Return parameter character containing TYPE NAME
  !!
  elemental function type(self)
    class(xCylinder), intent(in) :: self
    character(nameLen)           :: type

    type = TYPE_NAME

  end function type

  !!
  !! Return string with definition of this surface
  !!
  pure subroutine getDef(self,string)
    class(xCylinder), intent(in)           :: self
    character(:),allocatable,intent(inout) :: string

    string = printSurfDef(TYPE_NAME, [self % radius, self % origin])

  end subroutine getDef

  !!
  !! Returns an axis alligned bouding box of surface -ve halfspace
  !!
  pure subroutine boundingBox(self,origin, halfwidth)
    class(xCylinder), intent(in)            :: self
    real(defReal), dimension(3),intent(out) :: origin
    real(defReal), dimension(3),intent(out) :: halfwidth

    origin       = self % origin
    halfwidth    = self % radius
    halfwidth(1) = INFINITY

  end subroutine boundingBox

  !!
  !! Override base type function to returns .false.
  !!
  function cannotBeBoundary(self) result(itCant)
    class(xCylinder), intent(in) :: self
    logical(defBool)       :: itCant

    itCant = .not.ACCEPTS_BC

  end function cannotBeBoundary

  !!
  !! Calculate distance to cylinder along direction u
  !!
  elemental subroutine distance(self, dist, idx, r, u)
    class(xCylinder), intent(in)    :: self
    real(defReal), intent(out)      :: dist
    integer(shortInt), intent(out)  :: idx
    type(vector), intent(in)        :: r
    type(vector), intent(in)        :: u
    real(defReal)                   :: yBar, zBar
    real(defReal)                   :: k, c, a
    real(defReal)                   :: discriminant

    ! Set index
    idx = self % myIdx()

    yBar = r % v(2) - self % origin(2)
    zBar = r % v(3) - self % origin(3)

    k = yBar * u % v(2) + zBar * u % v(3)
    a = ONE - u % v(1) * u % v(1)
    c = yBar*yBar + zBar*zBar - self % rSquared
    discriminant = k*k - a*c

    ! No intersection
    if ((a == ZERO).OR.(discriminant < ZERO)) then
      dist = INFINITY
      return
    ! Particle is on the cylinder - the distance will be either positive or negative
    else if (abs(c) < surface_tol) then
      if (k >= ZERO) then
        dist = INFINITY
        return
      else
        dist = (-k + sqrt(discriminant))/a
        return
      end if
    ! Particle is inside - take solution with + before sqrt
    else if (c < ZERO) then
      dist = (-k + sqrt(discriminant))/a
      return
    ! Particle is outside - both distances are either positive or negative
    else
      dist = (-k -sqrt(discriminant))/a
      if (dist < ZERO) dist = INFINITY
      return
    end if

  end subroutine distance

  !!
  !! Return normal to the cylinder
  !!
  elemental function normalVector(self, r) result(normal)
    class(xCylinder), intent(in)   :: self
    type(vector), intent(in)       :: r
    type(vector)                   :: normal

    normal = [ZERO, TWO*(r % v(2) - self%origin(2)), TWO*(r % v(3) - self%origin(3))]

  end function normalVector

  !!
  !! Set boundary conditions for an xCylinder: may only be vacuum
  !!
  subroutine setBoundaryConditions(self, BC)
    class(xCylinder), intent(inout)             :: self
    integer(shortInt), dimension(:), intent(in) :: BC
    character(100),parameter :: Here ='setBoundaryConditions (xCylinder_class.f90)'

    if (any(BC /= vacuum)) then
      call fatalError(Here,'Cylinder boundaries may only be vacuum')

    else
      self % BC = vacuum

    end if

  end subroutine setBoundaryConditions

  !!
  !! Apply boundary transformation
  !!
  subroutine boundaryTransform(self, r, u)
    class(xCylinder), intent(in)     :: self
    type(vector), intent(inout)      :: r
    type(vector), intent(inout)      :: u
    character(100),parameter         :: Here ='boundaryTransform (xCylinder_class.f90)'

    select case(self % BC)
      case(vacuum)
        return

      case(noBC)
        call fatalError(Here, 'This surface has no BC')

      case default
        call fatalError(Here, 'Unsuported BC')

      end select
  end subroutine boundaryTransform

    
end module xCylinder_class
