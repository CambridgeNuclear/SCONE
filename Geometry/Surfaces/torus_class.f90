module torus_class

  use numPrecision
  use universalVariables, only : SURF_TOL, INF, X_AXIS, Y_AXIS, Z_AXIS
  use genericProcedures,  only : fatalError, numToChar
  use dictionary_class,   only : dictionary
  use surface_inter,      only : kill_super => kill

  implicit none
  private

  !!
  !! Axis aligned elliptic torus
  !!
  !! F(r) = (ri - ri0)^2 / B^2 + (SQRT[(rj - rjo)^2 + (rk - rko)^2] - A)^2 / C^2 - 1
  !!
  !! Where i,j,k (i /= j /= k) can be any of x,y & z axis
  !! A is the major radius. B and C correspond to the minor radii which can
  !! generally be eccentric. ri0, rj0, and rk0 are the origin of the torus
  !!
  !! Three diffrent types are avaliable
  !!   xTorus -> aligned with X-axis
  !!   yTorus -> aligned with Y-axis
  !!   zTorus -> aligned with Z-axis
  !!
  !! Surface tolerance: 2 * R * SURF_TOL ???
  !!
  !! Sample dictionary input:
  !!   tor { type xTorus; // could be yTorus or zTorus as well
  !!         id 3;
  !!         origin (0.0 0.0 0.0);
  !!         coeffs (5 10 20; }  // correspond to A B C above
  !!
  !! Private Members:
  !!   axis   -> Index of an alignment axis in {X_AXIS, Y_AXIS, Z_AXIS}
  !!   plane  -> Indices of axes in plane of torus {X_AXIS, Y_AXIS, Z_AXIS}\{axis}
  !!   origin -> Origin of the torus
  !!   r      -> Toroidal radii A, B, and C
  !!   r_sq   -> Square of radius r (r^2)
  !!
  !! Interface:
  !!   surface interface
  !!
  type, public, extends(surface) :: torus
    private
    integer(shortInt)               :: axis   = 0
    integer(shortInt), dimension(2) :: plane  = [0, 0]
    real(defReal), dimension(3)     :: origin = ZERO
    real(defReal), dimension(3)     :: r      = ZERO
    real(defReal)                   :: r_sq   = ZERO
  contains
    ! Superclass procedures
    procedure :: myType
    procedure :: init
    procedure :: boundingBox
    procedure :: evaluate
    procedure :: distance
    procedure :: going
    procedure :: kill

    ! Local procedures
    procedure :: build

  end type torus

contains

  !!
  !! Return surface type name
  !!
  !! See surface_inter for more details
  !!
  pure function myType(self) result(str)
    class(torus), intent(in)  :: self
    character(:), allocatable :: str

    select case(self % axis)
      case(X_AXIS)
        str = 'xTorus'

      case(Y_AXIS)
        str = 'yTorus'

      case(Z_AXIS)
        str = 'zTorus'

      case default
        str = 'unknown torus'

    end select

  end function myType

  !!
  !! Initialise torus from a dictionary
  !!
  !! See surface_inter for more details
  !!
  !! Errors:
  !!   fatalError if radius or id < 0.
  !!
  subroutine init(self, dict)
    class(torus), intent(inout)              :: self
    class(dictionary), intent(in)            :: dict
    integer(shortInt)                        :: id
    real(defReal), dimension(:), allocatable :: origin, radii
    character(nameLen)                       :: name
    character(100), parameter :: Here = 'init (torus_class.f90)'

    ! Get from dictionary
    call dict % get(id, 'id')
    call dict % get(origin, 'origin')
    call dict % get(radii, 'coeffs')
    call dict % get(name, 'type')

    ! Ensure positive id
    if (id <= 0) call fatalError(Here,'Surface id must be positive')

    ! Check origin size
    if (size(origin) /= 3) then
      call fatalError(Here,'Origin needs to have size 3. Has: '//numToChar(size(origin)))
    end if
    
    ! Check radii size
    if (size(radii) /= 3) then
      call fatalError(Here,'coeffs needs to have size 3. Has: '//numToChar(size(radii)))
    end if

    ! Ensure radii positivity
    if (any(radii < 0)) call fatalError(Here,'Torus coeffs must be positive')

    ! Other necessary checks on radii?

    self % r = radii
    self % origin = origin
    call self % setID(id)

    ! Build torus
    call self % build(name)

    ! Set surface tolerance
    call self % setTol( TWO * self % r(1) * SURF_TOL)

  end subroutine init

  !!
  !! Build torus from components
  !!
  !! Args:
  !!   id [in]   -> Surface ID
  !!   type [in] -> Torus type {'xTorus', 'yTorus' or 'zTorus'}
  !!
  !! Errors:
  !!   fatalError if invalid torus type
  !!
  subroutine build(self, type)
    class(torus), intent(inout) :: self
    character(*), intent(in)    :: type
    character(100), parameter :: Here = 'build (torus_class.f90)'

    ! Select type of torus
    select case(type)
      case('xTorus')
        self % axis = X_AXIS
        self % plane = [Y_AXIS, Z_AXIS]

      case('yTorus')
        self % axis = Y_AXIS
        self % plane = [X_AXIS, Z_AXIS]

      case('zTorus')
        self % axis = Z_AXIS
        self % plane = [X_AXIS, Y_AXIS]

      case default
        call fatalError(Here, 'Unknown type of torus: '//type)

    end select

  end subroutine build

  !!
  !! Return axis-aligned bounding box for the surface
  !!
  !! See surface_inter for details
  !!
  pure function boundingBox(self) result(aabb)
    class(torus), intent(in)    :: self
    real(defReal), dimension(6) :: aabb

    ! Top and bottom in a plane
    aabb(self % plane) = self % origin(self % plane) - [self % r, self % r]
    aabb(3 + self % plane) = self % origin(self % plane) + [self % r, self % r]

    ! Axis
    aabb(self % axis) = -INF
    aabb(3 + self % axis) = INF

  end function boundingBox

  !!
  !! Evaluate surface expression c = F(r)
  !!
  !! See surface_inter for details
  !!
  pure function evaluate(self, r) result(c)
    class(torus), intent(in)                :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: c

    associate (x => r(self % axis), x0 => self % origin(self % axis), &
                    y => r(self % plane(1)), y0 => self % origin(self % plane(1)), &
                    z => r(self % plane(2)), z0 => self % origin(self % plane(2)))
      c = (x - x0)**2 / self % r2(1) + &
              (sqrt((y - y0)**2 + (z - z0)**2) - self % r(1)) / self % r(3) - 1

    end associate

  end function evaluate

  !!
  !! Return distance to the surface
  !!
  !! See surface_inter for details
  !!
  !! Solves quadratic intersection equation
  !!   ad^2 + 2kd + c = 0
  !!   c = F(r)
  !!   k = (r1-x0)u1 + (r2-y0)u2
  !!   a = u1^2 + u2^2
  !!
  pure function distance(self, r, u) result(d)
    class(torus), intent(in)                :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    real(defReal)                           :: d, a
    real(defReal)                           :: c, k, delta

    ! Calculate quadratic components in the plane
    c = self % evaluate(r)
    k = dot_product(r(self % plane) - self % origin(self % plane) , u(self % plane))
    a = ONE - u(self % axis)**2
    delta = k*k - a*c  ! Technically delta/4

    ! Calculate the distance
    if (delta < ZERO .or. a == ZERO) then ! No intersection
      d = INF

    else if (abs(c) < self % surfTol()) then ! Point at a surface
      if ( k >= ZERO) then
        d = INF
      else
        d = -k + sqrt(delta)
        d = d/a
      end if

    else if (c < ZERO) then ! Point inside the surface
      d = -k + sqrt(delta)
      d = d/a

    else ! Point outside the surface
      d = -k - sqrt(delta)
      d = d/a
      if (d <= ZERO) d = INF

    end if

    ! Cap distance at Infinity
    d = min(d, INF)

  end function distance

  !!
  !! Returns TRUE if particle is going into +ve halfspace
  !!
  !! See surface_inter for details
  !!
  pure function going(self, r, u) result(halfspace)
    class(torus), intent(in)                :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    logical(defBool)                        :: halfspace
    real(defReal), dimension(2)             :: rp, up

    rp = r(self % plane) - self % origin(self % plane)
    up = u(self % plane)

    halfspace = dot_product(rp , up ) >= ZERO

  end function going

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(torus), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Local
    self % axis   = 0
    self % plane   = 0
    self % origin = ZERO
    self % r      = ZERO
    self % r_sq   = ZERO

  end subroutine kill

end module torus_class
