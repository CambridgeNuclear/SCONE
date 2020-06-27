module aPlane_class

  use numPrecision
  use universalVariables, only : X_AXIS, Y_AXIS, Z_AXIS, INF
  use genericProcedures,  only : fatalError
  use dictionary_class,   only : dictionary
  use quadSurface_inter,  only : quadSurface
  use surface_inter,      only : kill_super => kill
  implicit none
  private

  !!
  !! (Axis) Plane
  !!
  !! Planes with a normal which is one of the principal axis (x,y or z)
  !!
  !! F(r) = (r1 -a0) = c
  !!
  !! Surface tolerance: SURF_TOL
  !!
  !! Sample dictionary input:
  !!  x { type xPlane; id 16; x0 3.0;}
  !!  y { type yPlane; id 19, y0 -1.2;}
  !!
  !! Private members:
  !!   axis -> Axis indentifier
  !!   a0   -> Offset. Position of the plane on the axis
  !!
  !! Interface:
  !!   surface interface
  !!
  type, public, extends(quadSurface) :: aPlane
    private
    integer(shortInt) :: axis = -7
    real(defReal)     :: a0   = ZERO
  contains
    ! Superclass procedures
    procedure :: myType
    procedure :: init
    procedure :: boundingBox
    procedure :: evaluate
    procedure :: distance
    procedure :: going
    procedure :: kill
  end type aPlane


contains

  !!
  !! Return surface type name
  !!
  !! See surface_inter for more details
  !!
  pure function myType(self) result(str)
    class(aPlane), intent(in)  :: self
    character(:), allocatable  :: str

    select case(self % axis)
      case(X_AXIS)
        str = 'xPlane'

      case(Y_AXIS)
        str = 'yPlane'

      case(Z_AXIS)
        str = 'zPlane'

      case default
        str = 'unknown aPlane'

    end select
  end function myType

  !!
  !! Initialise aPlane from a dictionary
  !!
  !! See surface_inter for more details
  !!
  !! Errors:
  !!   fatalError if id < 0.
  !!
  subroutine init(self, dict)
    class(aPlane), intent(inout)  :: self
    class(dictionary), intent(in) :: dict
    integer(shortInt)             :: id
    character(nameLen)            :: type
    character(100), parameter :: Here = 'init (sphere_class.f90)'

    ! Get from dictionary
    call dict % get(id, 'id')
    call dict % get(type,'type')

    ! Check values
    if (id < 1) then
      call fatalError(Here,'Invalid surface id provided. ID must be > 1')
    end if

    ! Load data
    call self % setID(id)

    select case(type)
      case('xPlane')
        self % axis = X_AXIS
        call dict % get(self % a0, 'x0')

      case('yPlane')
        self % axis = Y_AXIS
        call dict % get(self % a0, 'y0')

      case('zPlane')
        self % axis = Z_AXIS
        call dict % get(self % a0, 'z0')

      case default
        call fatalError(Here, 'Unknown type of axis plane: '//type)
    end select

  end subroutine init

  !!
  !! Return axix-align bounding box for the surface
  !!
  !! See surface_inter for details
  !!
  pure function boundingBox(self) result(aabb)
    class(aPlane), intent(in)   :: self
    real(defReal), dimension(6) :: aabb

    aabb(1:3) = -INF
    aabb(4:6) = INF

    aabb([self % axis, self % axis + 3]) = self % a0

  end function boundingBox

  !!
  !! Evaluate surface expression c = F(r)
  !!
  !! See surface_inter for details
  !!
  pure function evaluate(self, r) result(c)
    class(aPlane), intent(in)               :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: c

    c = r(self % axis) - self % a0

  end function evaluate


  !!
  !! Return distance to the surface
  !!
  !! See surface_inter for details
  !!
  pure function distance(self, r, u) result(d)
    class(aPlane), intent(in)              :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    real(defReal)                           :: d
    real(defReal)                           :: ua, ra

    ra = self % a0 - r(self % axis)
    ua = u(self % axis)

    if(abs(ra) < self % surfTol()) then ! Within surface tolerance
      d = INF

    else if (ua /= ZERO) then ! Normal case
      d = ra/ua

    else  ! Parallel to the plane
      d = INF

    end if

    ! Cap the distance
    if ( d <= ZERO .or. d > INF ) d = INF

  end function distance

  !!
  !! Returns TRUE if particle is going into +ve halfspace
  !!
  !! See surface_inter for details
  !!
  pure function going(self, r, u) result(halfspace)
    class(aPlane), intent(in)              :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    logical(defBool)                        :: halfspace
    real(defReal)                           :: ua

    ua = u(self % axis)
    halfspace = ua > ZERO

    ! Special case of parallel direction
    ! Partilce stays in ist current halfspace
    if (ua == ZERO) then
      halfspace = (r(self % axis) - self % a0) >= ZERO
    end if

  end function going

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(aPlane), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Local
    self % axis = -7
    self % a0 = ZERO

  end subroutine kill


end module aPlane_class
