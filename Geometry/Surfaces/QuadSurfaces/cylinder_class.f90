module cylinder_class

  use numPrecision
  use universalVariables, only : SURF_TOL, INF, X_AXIS, Y_AXIS, Z_AXIS
  use genericProcedures,  only : fatalError, numToChar, areEqual
  use dictionary_class,   only : dictionary
  use quadSurface_inter,  only : quadSurface
  use surface_inter,      only : kill_super => kill

  implicit none
  private

  !!
  !! Axis aligned cylinder
  !!
  !! F(r) = (ri - oi)^2 + (rj - oj)^2 - R^2 = 0
  !!
  !! Where i,j (i /= j) can be any of x,y & z axis
  !!
  !! Three diffrent types are avaliable
  !!   xCylinder -> aligned with X-axis
  !!   yCylinder -> aligned with Y-axis
  !!   zCylinder -> aligned with Z-axis
  !!
  !! Surface tolerance: 2 * R * SURF_TOL
  !!
  !! Sample dictionary input:
  !!   cyl { type xCylinder; // could be yCylinder or zCylinder as well
  !!         id 3;
  !!         origin (0.0 0.0 0.0);
  !!         radius 7.34; }
  !!
  !! Private Members:
  !!   axis   -> Index of an alignment axis in {X_AXIS, Y_AXIS, Z_AXIS}
  !!   plane  -> Indexes of axis in plane of cylinder {X_AXIS, Y_AXIS, Z_AXIS}\{axis}
  !!   origin -> Location of the middle of the cylinder (component in axis has no significance)
  !!   r      -> Cylinder radius
  !!   r_sq   -> Square of radius r (r^2)
  !!
  !! Interface:
  !!   surface interface
  !!
  type, public, extends(quadSurface) :: cylinder
    private
    integer(shortInt)               :: axis   = 0
    integer(shortInt), dimension(2) :: plane  = 0
    real(defReal), dimension(3)     :: origin = ZERO
    real(defReal)                   :: r      = ZERO
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

  end type cylinder

contains

  !!
  !! Return surface type name
  !!
  !! See surface_inter for more details
  !!
  pure function myType(self) result(str)
    class(cylinder), intent(in) :: self
    character(:), allocatable   :: str

    select case(self % axis)
      case(X_AXIS)
        str = 'xCylinder'

      case(Y_AXIS)
        str = 'yCylinder'

      case(Z_AXIS)
        str = 'zCylinder'

      case default
        str = 'unknown cylinder'

    end select

  end function myType

  !!
  !! Initialise cylinder from a dictionary
  !!
  !! See surface_inter for more details
  !!
  !! Errors:
  !!   fatalError if radius or id < 0.
  !!
  subroutine init(self, dict)
    class(cylinder), intent(inout)  :: self
    class(dictionary), intent(in)   :: dict
    integer(shortInt)                        :: id
    real(defReal), dimension(:), allocatable :: origin
    character(nameLen)                       :: name
    real(defReal)                            :: r
    character(100), parameter :: Here = 'init (cylinder_class.f90)'

    ! Get from dictionary
    call dict % get(id, 'id')
    call dict % get(r, 'radius')
    call dict % get(origin, 'origin')
    call dict % get(name, 'type')

    ! Check origin size
    if (size(origin) /= 3) then
      call fatalError(Here,'Origin needs to have size 3. Has: '//numToChar(size(origin)))
    end if

    ! Build cylinder
    call self % build(id, name, origin, r )

  end subroutine init

  !!
  !! Build cylinder from components
  !!
  !! Args:
  !!   id [in]   -> Surface ID
  !!   type [in] -> Cylinder type {'xCylinder', 'yCylinder' or 'zCylinder'}
  !!   origin [in] -> Cylinder origin
  !!   radius [in] -> Cylinder radius
  !!
  !! Errors:
  !!   fatalError if id or radius are -ve
  !!
  subroutine build(self, id, type, origin, radius)
    class(cylinder), intent(inout)          :: self
    integer(shortInt), intent(in)           :: id
    character(*), intent(in)                :: type
    real(defReal), dimension(3), intent(in) :: origin
    real(defReal), intent(in)               :: radius
    character(100), parameter :: Here = 'build (cylinder_class.f90)'

    ! Check values
   if (id < 1) then
     call fatalError(Here,'Invalid surface id provided. ID must be > 1')

   else if ( radius <= ZERO) then
     call fatalError(Here, 'Radius of cylinder must be +ve. Is: '//numToChar(radius))

   end if

   ! Select type of cylinder
   select case(type)
     case('xCylinder')
       self % axis = X_AXIS
       self % plane = [Y_AXIS, Z_AXIS]

     case('yCylinder')
       self % axis = Y_AXIS
       self % plane = [X_AXIS, Z_AXIS]

     case('zCylinder')
       self % axis = Z_AXIS
       self % plane = [X_AXIS, Y_AXIS]

     case default
       call fatalError(Here, 'Unknown type of cylinder: '//type)

   end select

   ! Load data
   self % r = radius
   self % r_sq = radius * radius
   self % origin = origin
   call self % setID(id)

   ! Set surface tolerance
   call self % setTol( TWO * self % r * SURF_TOL)

  end subroutine build

  !!
  !! Return axis-aligned bounding box for the surface
  !!
  !! See surface_inter for details
  !!
  pure function boundingBox(self) result(aabb)
    class(cylinder), intent(in)   :: self
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
    class(cylinder), intent(in)             :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: c
    real(defReal), dimension(2)             :: diff

    diff = r(self % plane) - self % origin(self % plane)
    c = dot_product(diff, diff) - self % r_sq

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
    class(cylinder), intent(in)             :: self
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
    if (delta < ZERO .or. areEqual(a, ZERO)) then ! No intersection
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
    class(cylinder), intent(in)              :: self
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
    class(cylinder), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Local
    self % axis   = 0
    self % plane   = 0
    self % origin = ZERO
    self % r      = ZERO
    self % r_sq   = ZERO

  end subroutine kill

end module cylinder_class
