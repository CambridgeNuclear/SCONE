module cone_class

  use numPrecision
  use universalVariables, only : SURF_TOL, INF, X_AXIS, Y_AXIS, Z_AXIS
  use genericProcedures,  only : fatalError, numToChar, areEqual
  use dictionary_class,   only : dictionary
  use quadSurface_inter,  only : quadSurface
  use surface_inter,      only : kill_super => kill

  implicit none
  private

  !!
  !! Axis aligned cone, truncated at hMin and hMax
  !!
  !! Surface description:
  !! F(r) = max{ (ri - vi)^2 + (rj - vj)^2 - t^2 * (rk - vk)^2;
  !!             (rk - vk) - hMax;
  !!            -(rk - vk) + hMin; }
  !!
  !! Where i,j,k are x,y & z axis, and the cone is aligned along the k axis; t is
  !! the tangent of the cone opening angle (defined as the angle between the axis
  !! and the cone surface). The point V(vx, vy, vz) is the cone vertex.
  !!
  !! The cone is truncated by two faces given by their coordinates along the axis
  !! of the cone (hMin and hMax) with the vertex being 0. The sign of hMin and hMax
  !! determines the orientation of the cone.
  !! NOTE: Both entries must have the same sign (single truncated cone), and
  !! hMin < hMax must be true.
  !!
  !! The opening angle of the cone must be provided in degrees, and be in the
  !! range 0 - 90 degrees (extremes excluded)
  !!
  !! Three different types are available
  !!   xCone -> aligned with X-axis
  !!   yCone -> aligned with Y-axis
  !!   zCone -> aligned with Z-axis
  !!
  !! Surface tolerance: SURF_TOL * meanRadius. meanRadius is the radius of the cone at
  !! mid height, between hMin and hMax. Ideally, the tolerance would be space dependent
  !! and vary with height.
  !!
  !! Sample dictionary input:
  !!   cone { type xCone; // could be yCone or zCone as well
  !!          id 3;
  !!          vertex (0.0 0.0 0.0);
  !!          angle 45;
  !!          hMin 0.0;
  !!          hMax 10.0; }
  !!
  !! Private Members:
  !!   axis      -> Index of an alignment axis in {X_AXIS, Y_AXIS, Z_AXIS}
  !!   plane     -> Indexes of axis in plane of cone {X_AXIS, Y_AXIS, Z_AXIS}\{axis}
  !!   vertex    -> Location of the vertex of the cone
  !!   tanSquare -> Square of the tangent of the opening angle
  !!   hMin      -> Cone lower boundary
  !!   hMax      -> Cone upper boundary
  !!
  !! Interface:
  !!   surface interface
  !!
  type, public, extends(quadSurface) :: cone
    private
    integer(shortInt)               :: axis   = 0
    integer(shortInt), dimension(2) :: plane  = 0
    real(defReal), dimension(3)     :: vertex = ZERO
    real(defReal)                   :: tanSquare = ZERO
    real(defReal)                   :: hMin      = ZERO
    real(defReal)                   :: hMax      = ZERO
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
  subroutine init(self, dict)
    class(cone), intent(inout)               :: self
    class(dictionary), intent(in)            :: dict
    integer(shortInt)                        :: id
    real(defReal), dimension(:), allocatable :: vertex
    character(nameLen)                       :: name
    real(defReal)                            :: angle, hMin, hMax
    character(100), parameter :: Here = 'init (cone_class.f90)'

    ! Get from dictionary
    call dict % get(id, 'id')
    call dict % get(name, 'type')
    call dict % get(angle, 'angle')
    call dict % get(vertex, 'vertex')
    call dict % get(hMin, 'hMin')
    call dict % get(hMax, 'hMax')

    ! Check values
    if (id < 1) call fatalError(Here,'Invalid surface id provided. ID must be > 0')

    if (size(vertex) /= 3) then
      call fatalError(Here,'Vertex needs to have size 3. It has size: '//numToChar(size(vertex)))
    end if

    if (angle <= ZERO .or. angle >= 90.0_defReal) then
      call fatalError(Here, 'Opening angle of cone must be in the range 0-90 degrees '//&
                            & '(extremes excluded). It is: '//numToChar(angle))
    end if

    if (hMin >= hMax) call fatalError(Here, 'hMin is greater than or equal to hMax.')
    if ((sign(hMin,hMax) /= hMin) .and. hMin /= ZERO .and. hMax /= ZERO) then
      call fatalError(Here, 'hMin and hMax have different signs.')
    end if

    ! Load properties
    self % vertex = vertex
    self % hMin = hMin
    self % hMax = hMax

    ! Build cone
    call self % build(id, name, angle)

  end subroutine init

  !!
  !! Build cone from components
  !!
  !! Args:
  !!   id [in]    -> Surface ID
  !!   type [in]  -> Cone type {'xCone', 'yCone' or 'zCone'}
  !!   angle [in] -> Cone opening angle
  !!
  !! Errors:
  !!   fatalError if the cone type provided is not recognised
  !!
  subroutine build(self, id, type, angle)
    class(cone), intent(inout)              :: self
    integer(shortInt), intent(in)           :: id
    character(*), intent(in)                :: type
    real(defReal), intent(in)               :: angle
    real(defReal)                           :: tangent, meanRadius
    character(100), parameter :: Here = 'build (cylinder_class.f90)'

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
   call self % setID(id)

   tangent = tan(angle * PI / 180.0_defReal)
   self % tanSquare = tangent * tangent

   ! Set surface tolerance depending on the cone radius at mid height
   meanRadius = (abs(self % hMin) + abs(self % hMax)) * HALF * tangent
   call self % setTol(SURF_TOL * meanRadius)

  end subroutine build

  !!
  !! Return axis-aligned bounding box for the surface
  !!
  !! See surface_inter for details
  !!
  pure function boundingBox(self) result(aabb)
    class(cone), intent(in)     :: self
    real(defReal), dimension(6) :: aabb
    real(defReal)               :: maxRadius

    ! Along the axis
    aabb(self % axis) = self % vertex(self % axis) + self % hMin
    aabb(self % axis + 3) = self % vertex(self % axis) + self % hMax

    ! On the plane of the bases
    maxRadius = max(abs(self % hMin), abs(self % hMax)) * sqrt(self % tanSquare)
    aabb(self % plane)     = self % vertex(self % plane) - maxRadius
    aabb(self % plane + 3) = self % vertex(self % plane) + maxRadius

  end function boundingBox

  !!
  !! Evaluate surface expression c = F(r)
  !!
  !! See surface_inter for details
  !!
  pure function evaluate(self, r) result(c)
    class(cone), intent(in)                 :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: c, cMin, cMax, cBase
    real(defReal), dimension(3)             :: diff

    ! Vector for the difference between position provided and vertex
    diff = r - self % vertex

    ! Evaluate the expression for the slanted surface of the cone
    c = dot_product(diff(self % plane), diff(self % plane)) - &
        self % tanSquare * diff(self % axis) ** 2

    ! Evaluate the expressions for the two bases of the cone and find the maximum
    cMin = - diff(self % axis) + self % hMin
    cMax = diff(self % axis) - self % hMax
    cBase = max(cMin, cMax)

    ! Find the overall maximum
    c = max(c, cBase)

  end function evaluate

  !!
  !! Return distance to the surface
  !!
  !! See surface_inter for details
  !!
  !! Solves quadratic intersection equation
  !!   a*d^2 + 2b*d + c = 0
  !!   c = cone surface expression
  !!   b = (r1 - v1)u1 + (r2 - v2)u2 - t^2(r3 - v3)u3
  !!   a = u1^2 + u2^2 - t^2 * u3^2
  !!
  pure function distance(self, r, u) result(dist)
    class(cone), intent(in)                 :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    real(defReal)                           :: dist
    real(defReal), dimension(3)             :: diff, orientation
    real(defReal)                           :: a, b, c, delta, d1, d2, near, far, &
                                               test_near, test_far
    logical(defBool)                        :: check
    real(defReal), parameter :: FP_MISS_TOL = ONE + 10.0_defReal * epsilon(ONE)

    ! Vector for the difference between position provided and vertex
    diff = r - self % vertex

    ! Explicitly specify cone orientation
    orientation(self % plane) = ZERO
    orientation(self % axis)  = sign(ONE, self % hMin)

    ! Calculate quadratic components in the cone
    a = dot_product(u(self % plane), u(self % plane)) - &
        self % tanSquare * u(self % axis)**2
    b = dot_product(diff(self % plane) , u(self % plane)) - &
        self % tanSquare * diff(self % axis) * u(self % axis)
    c = dot_product(diff(self % plane), diff(self % plane)) - &
        self % tanSquare * diff(self % axis) ** 2

    ! Calculate delta/4 (discriminant of the quadratic equation: a*d^2 + 2b*d + c = 0)
    delta = b * b - a * c

    ! Calculate the distances from the slanted surface
    if (abs(a) < epsilon(ONE)) then   ! One intersection: particle direction tangent to cone opening
      d1 = - HALF * c / b
      d2 = sign(INF, dot_product(u, orientation))

      ! Check that the intersection found is in the right hemicone. Otherwise set both to INF
      ! The sign of inf depends on the direction of the particle and the orientation of the cone
      check = dot_product((diff + d1 * u), orientation) > ZERO
      if (.not. check) d1 = d2

    elseif (delta < epsilon(ONE)) then        ! No intersection (should never happen)
      d1 = -INF
      d2 = INF

    else                              ! Two intersections
      d1 = - (b + sqrt(delta)) / a
      d2 = - (b - sqrt(delta)) / a

      ! Check that the intersections found are in the right hemicone. Otherwise set to INF
      ! The sign of inf depends on the direction of the particle and the orientation of the cone
      check = dot_product((diff + d1 * u), orientation) > ZERO
      if (.not. check) d1 = sign(INF, dot_product(u, orientation))

      check = dot_product((diff + d2 * u), orientation) > ZERO
      if (.not. check) d2 = sign(INF, dot_product(u, orientation))

    end if

    ! Save minimum and maximum distance from cone surface
    near = min(d1, d2)
    far  = max(d1, d2)

    ! Calculate the distances from the bases of the cone
    if (abs(u(self % axis)) > epsilon(ONE)) then    ! Normal intersection

      d1 = (self % hMax - diff(self % axis)) / u(self % axis)
      d2 = (self % hMin - diff(self % axis)) / u(self % axis)

    else                                            ! Particle parallel to axis: check location

      if (diff(self % axis) > self % hMin .and. diff(self % axis) < self % hMax) then
        ! Inside the cone: base intersection segment lies between -INF:+INF
        d1 = -INF
        d2 = INF
      else
        ! Outside the cone: base intersection segment doesn't exist
        d1 = -INF
        d2 = -INF
      end if

    end if

    ! Save minimum and maximum distance from cone bases
    test_near = min(d1, d2)
    test_far  = max(d1, d2)

    ! Get intersection between sets of distances
    far  = min(far, test_far)
    near = max(near, test_near)

    ! Choose correct distance
    if (far <= near * FP_MISS_TOL) then      ! There is no intersection
      dist = INF

    else if (abs(self % evaluate(r)) < self % surfTol()) then ! Point at the surface
      ! Choose distance with largest absolute value
      if (abs(far) >= abs(near)) then
        dist = far
      else
        dist = near
      end if

    else                                     ! Normal hit. Closest distance
      if (near <= ZERO) then
        dist = far
      else
        dist = near
      end if

    end if

    ! Cap the distance
    if (dist <= ZERO .or. dist > INF) then
      dist = INF
    end if

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
    real(defReal), dimension(3)             :: diff, norm
    real(defReal)                           :: proj

    diff = r - self % vertex

    ! Check the location of the particle, i.e., base or cone surface, to calculate
    ! the normal
    if (abs(diff(self % axis) - self % hMin) < self % surfTol()) then
      norm(self % axis)  = -ONE
      norm(self % plane) = ZERO

    elseif (abs(diff(self % axis) - self % hMax) < self % surfTol()) then
      norm(self % axis)  = ONE
      norm(self % plane) = ZERO

    else
      norm(self % plane) = diff(self % plane)
      norm(self % axis)  = - self % tanSquare * diff(self % axis)

    end if

    norm = norm/norm2(norm)
    proj = dot_product(norm,u)

    ! Parallel direction. Need to use position to determine halfspace.
    if (areEqual(proj, ZERO)) then
      halfspace = self % evaluate(r) >= ZERO
      return
    
    end if

    ! Determine halfspace
    halfspace = proj > ZERO

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
    self % tanSquare = ZERO
    self % hMin      = ZERO
    self % hMax      = ZERO

  end subroutine kill

end module cone_class
