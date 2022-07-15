module cone_class

  use numPrecision
  use universalVariables, only : SURF_TOL, INF, X_AXIS, Y_AXIS, Z_AXIS
  use genericProcedures,  only : fatalError, numToChar, dotProduct
  use dictionary_class,   only : dictionary
  use quadSurface_inter,  only : quadSurface
  use surface_inter,      only : kill_super => kill

  implicit none
  private

  !!
  !! Axis aligned cone
  !!
  !! F(r) = (ri - oi)^2 + (rj - oj)^2 - t^2 * (rk - ok)^2 = 0
  !!
  !! Where i,j,k are x,y & z axis, and t is the tangent of the cone opening angle
  !! The cone is aligned along the k axis
  !!
  !! Three diffrent types are avaliable
  !!   xCone -> aligned with X-axis
  !!   yCone -> aligned with Y-axis
  !!   zCone -> aligned with Z-axis
  !!
  !! Surface tolerance: 2 * R * SURF_TOL
  !!
  !! Sample dictionary input:
  !!   cone { type xCone; // could be yCone or zCone as well
  !!         id 3;
  !!         vertex (0.0 0.0 0.0);
  !!         tangent 0.85;
  !!         orientation 1; }
  !!
  !! Private Members:
  !!   axis   -> Index of an alignment axis in {X_AXIS, Y_AXIS, Z_AXIS}
  !!   plane  -> Indexes of axis in plane of cone {X_AXIS, Y_AXIS, Z_AXIS}\{axis}
  !!   vertex -> Location of the vertex of the cone
  !!   t      -> Tangent of the opening angle
  !!
  !! Interface:
  !!   surface interface
  !!
  type, public, extends(quadSurface) :: cone
    private
    integer(shortInt)               :: axis   = 0
    integer(shortInt)               :: dir    = 0
    integer(shortInt), dimension(2) :: plane  = 0
    real(defReal), dimension(3)     :: vertex = ZERO
    real(defReal)                   :: t_sq   = ZERO
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
  end type cone

contains

  !!
  !! Return surface type name
  !!
  !! See surface_inter for more details
  !!
  pure function myType(self) result(str)
    class(cone), intent(in)   :: self
    character(:), allocatable :: str

    select case(self % axis)
      case(X_AXIS)
        str = 'xCone'

      case(Y_AXIS)
        str = 'yCone'

      case(Z_AXIS)
        str = 'zCone'

      case default
        str = 'unknown cone'

    end select

  end function myType

  !!
  !! Initialise cone from a dictionary
  !!
  !! See surface_inter for more details
  !!
  !! Errors:
  !!   fatalError if radius or id < 0.
  !!
  subroutine init(self, dict)
    class(cone), intent(inout)      :: self
    class(dictionary), intent(in)   :: dict
    integer(shortInt)                        :: id, dir
    real(defReal), dimension(:), allocatable :: vertex
    character(nameLen)                       :: name
    real(defReal)                            :: t
    character(100), parameter :: Here = 'init (cone_class.f90)'

    ! Get from dictionary
    call dict % get(id, 'id')
    call dict % get(t, 'tangent')
    call dict % get(vertex, 'vertex')
    call dict % get(name, 'type')
    call dict % get(dir, 'orientation')

    ! Check origin size
    if (size(vertex) /= 3) then
      call fatalError(Here,'Vertex needs to have size 3. Has: '//numToChar(size(vertex)))
    end if

    ! Build cylinder
    call self % build(id, name, vertex, t, dir)

  end subroutine init

  !!
  !! Build cone from components
  !!
  !! Args:
  !!   id [in]   -> Surface ID
  !!   type [in] -> Cone type {'xCone', 'yCone' or 'zCone'}
  !!   vertex [in]  -> Cone vertex
  !!   tangent [in] -> Cone opening angle tangent
  !!   dir [in]     -> Orientation of the cone
  !!
  !! Errors:
  !!   fatalError if id or radius are -ve
  !!
  subroutine build(self, id, type, vertex, tangent, dir)
    class(cone), intent(inout)              :: self
    integer(shortInt), intent(in)           :: id
    character(*), intent(in)                :: type
    real(defReal), dimension(3), intent(in) :: vertex
    real(defReal), intent(in)               :: tangent
    integer(shortInt), intent(in)           :: dir
    character(100), parameter :: Here = 'build (cylinder_class.f90)'

    ! Check values
   if (id < 1) then
     call fatalError(Here,'Invalid surface id provided. ID must be > 1')

   else if ( tangent <= ZERO) then
     call fatalError(Here, 'Tangent of cone must be +ve. It is: '//numToChar(tangent))

   else if (dir /= 1 .and. dir /= -1) then
     call fatalError(Here, 'Orientation of cone must be 1 or -1. It is: '//numToChar(dir))

   end if

   ! Select type of cylinder
   select case(type)
     case('xCone')
       self % axis = X_AXIS
       self % plane = [Y_AXIS, Z_AXIS]

     case('yCone')
       self % axis = Y_AXIS
       self % plane = [X_AXIS, Z_AXIS]

     case('zCone')
       self % axis = Z_AXIS
       self % plane = [X_AXIS, Y_AXIS]

     case default
       call fatalError(Here, 'Unknown type of cone: '//type)

   end select

   ! Load data
   self % t_sq   = tangent * tangent
   self % vertex = vertex
   self % dir = dir
   call self % setID(id)

   ! Set surface tolerance
   call self % setTol(SURF_TOL)

  end subroutine build

  !!
  !! Return axis-aligned bounding box for the surface
  !!
  !! See surface_inter for details
  !!
  pure function boundingBox(self) result(aabb)
    class(cone), intent(in)   :: self
    real(defReal), dimension(6) :: aabb

    ! Top and bottom in a plane
    aabb(self % plane) = -INF
    aabb(3 + self % plane) = INF

    ! Axis
    if (self % dir == 1) then
      aabb(self % axis) = self % vertex(self % axis)
      aabb(3 + self % axis) = INF
    else
      aabb(self % axis) = - INF
      aabb(3 + self % axis) = self % vertex(self % axis)
    end if

  end function boundingBox

  !!
  !! Evaluate surface expression c = F(r)
  !!
  !! See surface_inter for details
  !!
  pure function evaluate(self, r) result(c)
    class(cone), intent(in)                 :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: c
    real(defReal), dimension(3)             :: diff

    diff = r - self % vertex
    c = dot_product(diff(self % plane), diff(self % plane)) - self % t_sq * diff(self % axis) ** 2

  end function evaluate

  !!
  !! Return distance to the surface
  !!
  !! See surface_inter for details
  !!
  !! Solves quadratic intersection equation
  !!   ad^2 + 2kd + c = 0
  !!   c = F(r)
  !!   k = (r1-x0)u1 + (r2-y0)u2 - t^2(r3-z0)u3
  !!   a = u1^2 + u2^2 - t^2 * u3^2
  !!
  pure function distance(self, r, u) result(d)
    class(cone), intent(in)             :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    real(defReal)                           :: d, a, c, k, delta, b

    ! Calculate quadratic components in the plane
    c = self % evaluate(r)
    k = dot_product(r(self % plane) - self % vertex(self % plane) , u(self % plane)) - &
        self % t_sq * r(self % axis) - self % vertex(self % axis)
    a = dot_product(u(self % plane), u(self % plane)) - self % t_sq * u(self % axis)**2
    delta = k*k - a*c  ! Technically delta/4

    ! Calculate the distance
    if (delta < ZERO .or. a == ZERO) then ! No intersection
      d = INF

    else if (abs(c) < self % surfTol()) then ! Point at a surface
      if ( k >= ZERO) then
        d = d/a
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
    class(cone), intent(in)                 :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    logical(defBool)                        :: halfspace
    real(defReal), dimension(2)             :: rp, up

  end function going

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(cone), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Local
    self % axis   = 0
    self % plane  = 0
    self % vertex = ZERO
    self % t_sq   = ZERO

  end subroutine kill

end module cone_class
