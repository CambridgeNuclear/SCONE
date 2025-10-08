module wedge_class

  use numPrecision
  use universalVariables
  use genericProcedures,  only : fatalError, numToChar, swap
  use dictionary_class,   only : dictionary
  use surface_inter,      only : surface, kill_super => kill

  implicit none
  private

  !!
  !! Axis Aligned wedge
  !!
  !! Wedge with two isosceles triangular bases, parallel between each other, and aligned
  !! with the x, y or z axis. The input type has to be ``xWedge``, ``yWedge`` or ``zWedge``.
  !! The wedge bases are characterised by a half opening angle (theta); the wedge can also
  !! be arbitrarily rotated around its axis by phi.
  !!
  !! F(r) = maxval( abs(r_ax - o_ax) - hw,
  !!                n1 \cdot (r_p - o_p),
  !!                n2 \cdot (r_p - o_p),
  !!                n3 \cdot (r_p - o_p) - h,
  !!              )
  !!
  !! The nomenclature for the 5 faces is:
  !!
  !!   - face1 and face2 are the two slanted faces defined by the opening angle:
  !!     face1 is the face rotated by -opening compared to the triangle altitude and
  !!     face2 is rotated by +opening.
  !!   - face3 refers to the face in front of the axis of the wedge
  !!   - -base is the lower base
  !!   - base is the upper base
  !!
  !! Where: hw -> halfwidth vector, o -> origin position
  !!        h  -> triangle height,  n  -> normal vectors
  !!        maxval -> maximum element
  !!
  !! Surface Tolerance: SURF_TOL
  !!
  !! Sample Dictionary Input:
  !!   aab { type zWedge; id 92; origin (0.0 0.0 9.0); halfwidth 5.0; altitude 2.0;
  !!         opening 30.0; rotation 0.0; }
  !!
  !! Boundary Conditions:
  !!   BC order: face1, face2, face3, -base, base
  !!   Each face can have diffrent BC. Any combination is supported with co-ordinate transform.
  !!
  !! Private Members:
  !!   origin -> Position of the middle of the wedge edge (or axis)
  !!   halfwidth -> Half-length of the wedge edge in the direction of the axis (must be > 0.0)
  !!   height -> Altitude of the triangular face of the wedge (must be > 0.0)
  !!   theta  -> Half-angle opening of the triangular face of the wedge (must be between
  !!             0.0 and 90.0 degrees excluded)
  !!   phi    -> Rotation angle around the edge of the wedge. The rotation angle is with respect
  !!             to the axis: +y for a xWedge; +x for a yWedge and zWedge. It must be between
  !!             0.0 and 360.0 degrees excluded
  !!   BC     -> Boundary conditions flags for each face (face1, face2, face3, -base, base)
  !!
  !! Interface:
  !!   surface interface
  !!
  type, public, extends(surface) :: wedge
    private
    real(defReal) :: halfwidth = ZERO
    real(defReal) :: height = ZERO
    real(defReal) :: theta  = ZERO
    real(defReal) :: phi    = ZERO
    real(defReal), dimension(2)     :: norm1   = ZERO
    real(defReal), dimension(2)     :: norm2   = ZERO
    real(defReal), dimension(2)     :: norm3   = ZERO
    real(defReal), dimension(3)     :: origin  = ZERO
    integer(shortInt)               :: axis    = 0
    integer(shortInt), dimension(2) :: plane   = 0
    integer(shortInt), dimension(5) :: BC      = VACUUM_BC

  contains
    ! Superclass procedures
    procedure :: myType
    procedure :: init
    procedure :: boundingBox
    procedure :: evaluate
    procedure :: distance
    procedure :: going
    procedure :: kill
    procedure :: setBC
    procedure :: explicitBC
    procedure :: transformBC

    ! Local procedures
    procedure :: build

  end type wedge


contains

  !!
  !! Return surface type name
  !!
  !! See surface_inter for more details
  !!
  pure function myType(self) result(str)
    class(wedge), intent(in) :: self
    character(:), allocatable  :: str

    select case(self % axis)
      case(X_AXIS)
        str = 'xWedge'

      case(Y_AXIS)
        str = 'yWedge'

      case(Z_AXIS)
        str = 'zWedge'

      case default
        str = 'unknown wedge'

    end select

  end function myType

  !!
  !! Initialise wedge from a dictionary
  !!
  !! See surface_inter for more details
  !!
  !! Errors:
  !!
  subroutine init(self, dict)
    class(wedge), intent(inout)     :: self
    class(dictionary), intent(in)   :: dict
    integer(shortInt)               :: id, N
    real(defReal)                   :: value, theta, phi
    character(nameLen)              :: type
    real(defReal), dimension(:), allocatable :: temp
    character(100), parameter :: Here = 'init (wedge_class.f90)'

    ! Load id
    call dict % get(id,'id')
    if (id <= 0) call fatalError(Here, 'ID must be <=0. Is: '//numToChar(id))
    call self % setID(id)

    ! Load origin
    call dict % get(temp,'origin')
    N = size(temp)
    if (N /= 3) call fatalError(Here,'origin must have size 3. Has: '//numToChar(N))
    self % origin = temp

    ! Load halfwidth
    call dict % get(value,'halfwidth')
    if (value <= ZERO) call fatalError(Here, 'halfwidth cannot be -ve.')
    self % halfwidth = value

    ! Load triangular face altitude
    call dict % get(value,'altitude')
    if (value <= ZERO) call fatalError(Here, 'altitude cannot be -ve.')
    self % height = value

    ! Load type - defines orientation
    call dict % get(type, 'type')

    ! Load triangle opening angle
    call dict % get(theta,'opening')

    ! Load rotation angle with respect to reference axis
    call dict % getOrDefault(phi,'rotation',ZERO)

    ! Build wedge
    call self % build(type, theta, phi)

  end subroutine init

  !!
  !! Build wedge from components
  !!
  !! Args:
  !!   type [in]  -> Wedge type {'xWedge', 'yWedge' or 'zWedge'}
  !!   theta [in] -> Wedge opening angle
  !!   phi [in]   -> Wedge rotation angle with respect to plane(1)
  !!
  !! Errors:
  !!   fatalError if the wedge type provided is not recognised
  !!
  subroutine build(self, type, theta, phi)
    class(wedge), intent(inout)         :: self
    character(*), intent(in)            :: type
    real(defReal), intent(in)           :: theta
    real(defReal), intent(in)           :: phi
    character(100), parameter :: Here = 'build (wedge_class.f90)'

   ! Select type of wedge
   select case(type)
   case('xWedge')
       self % axis = X_AXIS
       self % plane = [Y_AXIS, Z_AXIS]

     case('yWedge')
       self % axis = Y_AXIS
       self % plane = [X_AXIS, Z_AXIS]

     case('zWedge')
       self % axis = Z_AXIS
       self % plane = [X_AXIS, Y_AXIS]

     case default
       call fatalError(Here, 'Unknown type of wedge: '//type)

   end select

   ! Save opening angle
   if (theta <= ZERO .or. theta >= 90.0_defReal) then
     call fatalError(Here, 'Opening angle of wedge must be in the range 0-90 degrees '//&
                           & '(extremes excluded). It is: '//numToChar(theta))
   end if
   self % theta = theta * PI / 180.0_defReal

   ! Save rotation angle
   if (phi <= ZERO .or. phi >= 360.0_defReal) then
     call fatalError(Here, 'Rotation angle of wedge must be in the range 0-360 degrees '//&
                           & '(extremes excluded). It is: '//numToChar(phi))
   end if
   self % phi = phi * PI / 180.0_defReal

   ! Compute normals of the three triangle faces
   self % norm1 = [sin(self % phi - self % theta), -cos(self % phi - self % theta)]
   self % norm2 = [-sin(self % phi + self % theta), cos(self % phi + self % theta)]
   self % norm3 = [cos(self % phi), sin(self % phi)]

   ! Load data
   call self % setTol(SURF_TOL * self % halfwidth)

  end subroutine build

  !!
  !! Return axis-aligned bounding box for the surface
  !!
  !! See surface_inter for details
  !!
  pure function boundingBox(self) result(aabb)
    class(wedge), intent(in)    :: self
    real(defReal), dimension(6) :: aabb
    real(defReal), dimension(2) :: p0, p1, p2, tg1, tg2

    ! Along the axis
    aabb(self % axis)     = self % origin(self % axis) - self % halfwidth
    aabb(self % axis + 3) = self % origin(self % axis) + self % halfwidth

    ! Calculate triangle points
    tg1 = [-self % norm1(2), self % norm1(1)]
    tg2 = [self % norm2(2), -self % norm2(1)]
    p0 = self % origin(self % plane)
    p1 = p0 + self % height * tg1 / cos(self % theta)
    p2 = p0 + self % height * tg2 / cos(self % theta)

    ! On the 2D plane
    aabb(self % plane(1)) = min( p0(1), p1(1), p2(1) )
    aabb(self % plane(2)) = min( p0(2), p1(2), p2(2) )
    aabb(self % plane(1) + 3) = max( p0(1), p1(1), p2(1) )
    aabb(self % plane(2) + 3) = max( p0(2), p1(2), p2(2) )

  end function boundingBox

  !!
  !! Evaluate surface expression c = F(r)
  !!
  !! See surface_inter for details
  !!
  pure function evaluate(self, r) result(c)
    class(wedge), intent(in)                :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: c
    real(defReal), dimension(3)             :: diff
    real(defReal)                           :: cAxis, cBase, cSide1, cSide2

    ! Displacement from apex/origin
    diff = r - self % origin

    ! Axis slab function
    cAxis = abs(diff(self % axis)) - self % halfwidth

    ! Side functions
    cSide1 = dot_product(self % norm1, diff(self % plane))
    cSide2 = dot_product(self % norm2, diff(self % plane))

    ! Base function
    cBase = dot_product(self % norm3, diff(self % plane)) - self % height

    ! Overall implicit function: negative inside
    c = max(cAxis, cBase, cSide1, cSide2)

  end function evaluate

  !!
  !! Return distance to the surface
  !!
  !! See surface_inter for details
  !!
  !!
  pure function distance(self, r, u) result(dist)
    class(wedge), intent(in)                :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    real(defReal)                           :: dist
    real(defReal), dimension(3)             :: diff
    real(defReal)                           :: d1, d2, k, c, near, far, test_near, test_far
    real(defReal), parameter :: FP_MISS_TOL = ONE + 10.0_defReal * epsilon(ONE)

    ! Displacement from origin
    diff = r - self % origin

    ! Initialise intersection set
    near = -INF
    far  = INF

    ! Calculate distance from the three faces individually
    ! Face 1

    k = dot_product(u(self % plane), self % norm1)
    c = dot_product(diff(self % plane), self % norm1)

    if (abs(k) < epsilon(ONE)) then
      if (c > ZERO) then
        near = INF
        far  = -INF
      end if
    else
      d1 = -c/k
      if (k > ZERO) then
        far = min(far, d1)
      else
        near = max(near, d1)
      end if
    end if

    ! Face 2
    k = dot_product(u(self % plane), self % norm2)
    c = dot_product(diff(self % plane), self % norm2)

    if (abs(k) < epsilon(ONE)) then
      if (c > ZERO) then
        near = INF
        far  = -INF
      end if
    else
      d1 = -c/k
      if (k > ZERO) then
        far = min(far, d1)
      else
        near = max(near, d1)
      end if
    end if

    ! Face 3
    k = dot_product(u(self % plane), self % norm3)
    c = dot_product(diff(self % plane), self % norm3)

    if (abs(k) < epsilon(ONE)) then
      if (c > self % height) then
        near = INF
        far  = -INF
      end if
    else
      d1 = -(c - self % height)/k
      if (k > ZERO) then
        far = min(far, d1)
      else
        near = max(near, d1)
      end if
    end if

    ! Calculate the distances from the bases of the wedge
    if (abs(u(self % axis)) > epsilon(ONE)) then    ! Normal intersection

      d1 = (self % halfwidth - diff(self % axis)) / u(self % axis)
      d2 = (-self % halfwidth - diff(self % axis)) / u(self % axis)

      ! Save minimum and maximum distance from wedge bases
      test_near = min(d1, d2)
      test_far  = max(d1, d2)

    else  ! Particle parallel to axis: check location

      if (abs(diff(self % axis)) < self % halfwidth) then
        ! Inside the cone: base intersection segment lies between -INF:+INF
        test_near = -INF
        test_far = INF
      else
        ! Outside the cone: base intersection segment doesn't exist
        test_near = -INF
        test_far = -INF
      end if

    end if

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

    elseif (self % evaluate(r) < ZERO) then     ! Normal hit. Pick far if p is inside the wedge and vice versa
      dist = far

    else
      dist = near

    end if

    ! Cap the distance
    if (dist <= ZERO .or. dist > INF) dist = INF

  end function distance

  !!
  !! Returns TRUE if particle is going into +ve halfspace
  !!
  !! See surface_inter for details
  !!
  !! Works by:
  !!  1) Determine a plane in which direction is the closest
  !!  2) Use normal for this plane and project distance
  !!  3) Determinie halfspace based on sign of the projection
  !!
  !! Note:
  !!   For parallel direction halfspace is asigned by the sign of `evaluate` result.
  !!
  pure function going(self, r, u) result(halfspace)
    class(wedge), intent(in)                :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    logical(defBool)                        :: halfspace
    real(defReal), dimension(3)             :: diff, norm
    real(defReal)                           :: cSide1, cSide2, cBase, proj

    ! Displacement from origin
    diff = r - self % origin

    ! Helpers to evaluate wedge
    cSide1 = dot_product(self % norm1, diff(self % plane))
    cSide2 = dot_product(self % norm2, diff(self % plane))
    cBase = dot_product(self % norm3, diff(self % plane)) - self % height

    ! Initialise surface normal
    norm = ZERO

    ! Identify surface
    if (abs(diff(self % axis) - self % halfwidth) < self % surfTol()) then
      norm(self % axis)  = ONE
    elseif (abs(diff(self % axis) + self % halfwidth) < self % surfTol()) then
      norm(self % axis)  = -ONE
    elseif (abs(cSide1) < self % surfTol()) then
      norm(self % plane) = self % norm1
    elseif (abs(cSide2) < self % surfTol()) then
      norm(self % plane) = self % norm2
    else
      norm(self % plane) = self % norm3
    end if

    norm = norm/norm2(norm)
    proj = dot_product(norm, u)

    ! Determine halfspace
    halfspace = proj > ZERO

    ! Parallel direction
    ! Need to use position to determine halfspace
    if (abs(proj) < epsilon(ONE)) then
      halfspace = self % evaluate(r) >= ZERO
    end if

  end function going

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(wedge), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Local
    self % origin    = ZERO
    self % halfwidth = ZERO
    self % height    = ZERO
    self % theta     = ZERO
    self % phi       = ZERO
    self % BC    = VACUUM_BC
    self % norm1 = ZERO
    self % norm2 = ZERO
    self % norm3 = ZERO
    self % axis  = 0
    self % plane = 0

  end subroutine kill

  !!
  !! Set boundary conditions
  !!
  !! See surface_inter for details
  !!
  subroutine setBC(self, BC)
    class(wedge), intent(inout)                 :: self
    integer(shortInt), dimension(:), intent(in) :: BC
    integer(shortInt)                           :: i
    character(100),parameter :: Here = 'setBC (wedge_class.f90)'

    if (size(BC) /= 5) call fatalError(Here,'Wrong size of BC string. Must be 5')

    ! Load BC codes
    self % BC = BC

    ! Verify that all BC flags make sense
    do i = 1, 5
      select case(BC(i))
        case (VACUUM_BC)
          ! Do nothing, pass
        case (REFLECTIVE_BC, PERIODIC_BC)
          if (i == 3) call fatalError(Here,'Face opposite to wedge edge can only have vacuum BC')
        case default
          call fatalError(Here,'Unrecognised BC: '//numToChar(BC(i))//' in position: '//numToChar(i))
      end select
    end do

    ! Verify periodic BCs
    if (.not. all( (self % BC([1,4]) == PERIODIC_BC) .eqv. (self % BC([2,5]) == PERIODIC_BC) )) then
      call fatalError(Here,'Periodic BC need to be applied to oposite surfaces')
    end if

  end subroutine setBC

  !!
  !! Apply explicit BCs
  !!
  !! See surface_inter for details
  !!
  !! Note:
  !!   - Go through all directions in order to account for corners
  !!
  subroutine explicitBC(self, r, u)
    class(wedge), intent(in)                   :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u
    real(defReal), dimension(5,3)              :: norm
    real(defReal), dimension(3)                :: diff, normAxis
    real(defReal), dimension(5)                :: c
    real(defReal)                              :: angle, u1, u2, proj
    integer(shortInt)                          :: i, bc
    character(100), parameter :: Here = 'explicitBC (wedge_class.f90)'

    ! Displacement from origin
    diff = r - self % origin

    ! Helpers to evaluate wedge
    c    = ZERO
    c(1) = dot_product(self % norm1, diff(self % plane))
    c(2) = dot_product(self % norm2, diff(self % plane))
    c(3) = dot_product(self % norm3, diff(self % plane)) - self % height
    c(4) = diff(self % axis) + self % halfwidth
    c(5) = diff(self % axis) - self % halfwidth

    norm = ZERO
    norm(1, self % plane) = self % norm1
    norm(2, self % plane) = self % norm2
    norm(3, self % plane) = self % norm3
    norm(4, self % axis)  = -ONE
    norm(5, self % axis)  = ONE

    ! Loop over directions
    axis : do i  = 1, 5

      ! Skip if particle is well inside the domain
      if (abs(c(i)) > self % surfTol()) cycle axis

      ! Choose correct BC and normal
      bc = self % BC(i)
      normAxis = norm(i,:)

      ! Apply BC
      select case(bc)
        case (VACUUM_BC)
          ! Do nothing. Pass

        case (REFLECTIVE_BC)
          proj = dot_product(u, normAxis)
          u = u - TWO * proj * normAxis

        case (PERIODIC_BC)

          ! Recalculate displacement from origin
          diff = r - self % origin

          if (i == 1 .or. i == 2) then

            ! Pick the right rotation angle
            if (i == 1) then
              angle =  TWO * self % theta
            else
              angle = -TWO * self % theta
            end if

            r(self % plane(1)) = self % origin(self % plane(1)) + cos(angle) * diff(self % plane(1)) - &
                                 sin(angle) * diff(self % plane(2))
            r(self % plane(2)) = self % origin(self % plane(2)) + sin(angle) * diff(self % plane(1)) + &
                                 cos(angle) * diff(self % plane(2))
            u1 = u(self % plane(1))
            u2 = u(self % plane(2))
            u(self % plane(1)) = cos(angle) * u1 - sin(angle) * u2
            u(self % plane(2)) = sin(angle) * u1 + cos(angle) * u2

          elseif (i == 4 .or. i == 5) then
            r(self % axis) = r(self % axis) - TWO * sign(self % halfwidth, diff(self % axis))
          end if

        case default
          call fatalError(Here, 'Unrecognised BC: '// numToChar(bc))
      end select

    end do axis

  end subroutine explicitBC

  !!
  !! Apply co-ordinate transform BC
  !!
  !! See surface_inter for details
  !!
  !! Note:
  !!   - Order of transformations does not matter
  !!   - Calculate distance (in # of transformations) for each direction and apply them
  !!
  subroutine transformBC(self, r, u)
    class(wedge), intent(in)                   :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u
    real(defReal), dimension(2)                :: norm, diffPlane, proj
    integer(shortInt)                          :: transNumAx, t, bc, transNumRot
    real(defReal)                              :: a_bar, diff, diff1, diff2, angle, &
                                                  rotSign, u1, u2
    character(100), parameter :: Here = 'transformBC (wedge_class.f90)'

    ! Calculate halfwidth reduced by the surface_tolerance
    ! Necessary to capture particles at the boundary
    a_bar = self % halfwidth - self % surfTol()

    ! Calculate distance (in # of transformations) in each direction
    transNumAx = ceiling(abs(r(self % axis) - self % origin(self % axis)) / a_bar) / 2

    ! Loop over directions & number of transformations
    transAx : do t = 1, transNumAx
      ! Find position wrt origin
      diff = r(self % axis) - self % origin(self % axis)

      ! Choose correct BC
      if (diff < ZERO) then
        bc = self % BC(4)
      else
        bc = self % BC(5)
      end if

      ! Apply BC
      select case(bc)
        case (VACUUM_BC)
          ! Do nothing. Pass

        case (REFLECTIVE_BC)
          ! Calculate displacment and perform reflection
          r(self % axis) = r(self % axis) - TWO * (diff - sign(self % halfwidth, diff))
          u(self % axis) = -u(self % axis)

        case (PERIODIC_BC)
          ! Calculate displacement and perform translation
          r(self % axis) = r(self % axis) - TWO * sign(self % halfwidth, diff)

        case default
          call fatalError(Here, 'Unrecognised axial BC: '// numToChar(bc))
      end select

    end do transAx

    ! Rotations around the plane
    diffPlane = r(self % plane) - self % origin(self % plane)

    ! Compute angle relative to apex
    angle = atan2(diffPlane(2), diffPlane(1)) - self % phi
    if (abs(angle) > PI) angle = angle - sign(TWO_PI, angle)

    ! Compute how many full 2θ rotations we’re away from the principal wedge
    transNumRot = floor((angle + self % theta) / (TWO * self % theta))

    ! Loop over number of rotations
    rotLoop : do t = 1, abs(transNumRot)

      ! Recalculate relative position and angle as they might have changed
      diffPlane = r(self % plane) - self % origin(self % plane)
      angle = atan2(diffPlane(2), diffPlane(1)) - self % phi
      if (abs(angle) > PI) angle = angle - sign(TWO_PI, angle)

      ! Decide which boundary we crossed based on sign
      if (angle > 0) then
        bc = self % BC(2)   ! crossed face 2
        norm    = self % norm2
        rotSign = -ONE
      else
        bc = self % BC(1)   ! crossed face 1
        norm    = self % norm1
        rotSign = ONE
      end if

      select case (bc)
        case (VACUUM_BC)
          ! Do nothing (particle leaves domain)

        case (PERIODIC_BC)
          ! Rotate position and velocity back into domain
          angle = rotSign * TWO * self % theta

          ! Position
          diff1 = r(self % plane(1)) - self % origin(self % plane(1))
          diff2 = r(self % plane(2)) - self % origin(self % plane(2))
          r(self % plane(1)) = self % origin(self % plane(1)) + cos(angle) * diff1 - &
                               sin(angle) * diff2
          r(self % plane(2)) = self % origin(self % plane(2)) + sin(angle) * diff1 + &
                               cos(angle) * diff2

          ! Velocity
          u1 = u(self % plane(1))
          u2 = u(self % plane(2))
          u(self % plane(1)) = cos(angle) * u1 - sin(angle) * u2
          u(self % plane(2)) = sin(angle) * u1 + cos(angle) * u2

        case (REFLECTIVE_BC)
          ! Reflect position (push back inside domain)
          diffPlane = r(self % plane) - self % origin(self % plane)
          r(self % plane) = r(self % plane) - TWO * dot_product(diffPlane, norm) * norm

          ! Reflect velocity
          proj = dot_product(u(self % plane), norm)
          u(self % plane) = u(self % plane) - TWO * proj * norm

        case default
          call fatalError(Here, 'Unrecognised rotational BC: '//numToChar(bc))
      end select

    end do rotLoop

  end subroutine transformBC

end module wedge_class
