module wedge_class

  use numPrecision
  use universalVariables
  use genericProcedures,  only : fatalError, numToChar, swap
  use dictionary_class,   only : dictionary
  use compSurface_inter,  only : compSurface
  use surface_inter,      only : kill_super => kill

  implicit none
  private

  !!
  !! Axis Aligned wedge
  !!
  !! F(r) = maxval( abs(r_ax - o_ax) - hw,
  !!                n1 \cdot (r_p - o_p),
  !!                n2 \cdot (r_p - o_p),
  !!                n3 \cdot (r_p - o_p) - h,
  !!              )
  !!
  !! Where: hw -> halfwidth vector, o -> origin position
  !!        h  -> triangle height,  n  -> normal vectors
  !!        maxval -> maximum element
  !!
  !! Surface Tolerance: SURF_TOL
  !!
  !! Sample Dictionary Input:
  !!   aab { type zWedge; id 92; origin (0.0 0.0 9.0); halfwidth 5.0; height 2.0;
  !!         opening 30.0; rotation 0.0; }
  !!
  !! Boundary Conditions:
  !!   BC order: face1, face2, face3, ax_min, ax_max
  !!   Each face can have diffrent BC. Any combination is supported with co-ordinate transform.
  !!
  !! Private Members:
  !!   origin -> position of the middle of the wedge edge
  !!   halfwidth -> Halfwidth (half-length) of the wedge in the direction of the axis (must be > 0.0)
  !!   height -> Height of the triangular face of the wedge (must be > 0.0)
  !!   theta  -> Half angle opening of the triangular face of the wedge (must be between 0.0 and 90.0 degrees excluded)
  !!   phi    -> Rotation angle from the axis plane(1) to the height of the triangular face of the wedge
  !!             (must be between 0.0 and 360.0 degrees excluded)
  !!   BC     -> Boundary conditions flags for each face (ax_min, ax_max, face1, face2, face3)
  !!
  !! Interface:
  !!   surface interface
  !!
  type, public, extends(compSurface) :: wedge
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

    ! Load triangular face height
    call dict % get(value,'height')
    if (value <= ZERO) call fatalError(Here, 'height cannot be -ve.')
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
  function distance(self, r, u) result(dist)
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

    !print*, k, c
    !print*, '1 ', near, far

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

    !print*, k, c
    !print*, '2 ', near, far

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

    !print*, k, c
    !print*, '3 ', near, far

    ! Calculate the distances from the bases of the wedge
    if (abs(u(self % axis)) > epsilon(ONE)) then    ! Normal intersection

      d1 = (self % halfwidth - diff(self % axis)) / u(self % axis)
      d2 = (-self % halfwidth - diff(self % axis)) / u(self % axis)

      ! Save minimum and maximum distance from wedge bases
      test_near = min(d1, d2)
      test_far  = max(d1, d2)

    else                                            ! Particle parallel to axis: check location

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

    !print*, '4 ', test_near, test_far

    ! Get intersection between sets of distances
    far  = min(far, test_far)
    near = max(near, test_near)

    !print*, 'final ', near, far

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

    elseif (self % evaluate(r) < ZERO) then     ! Normal hit. Pick far if p is inside the wedge and viceversa
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
    class(wedge), intent(in)              :: self
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
    self % BC = VACUUM_BC

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
    class(wedge), intent(in)                     :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u
    integer(shortInt)                          :: ax, bc
    real(defReal)                              :: r0
    character(100), parameter :: Here = 'explicitBC (wedge_class.f90)'



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
    class(wedge), intent(in)                     :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u
    real(defReal), dimension(3)                :: a_bar
    integer(shortInt), dimension(3)            :: Ri
    integer(shortInt)                          :: ax, t, bc
    real(defReal)                              :: a0, d, r0
    character(100), parameter :: Here = 'transformBC (wedge_class.f90)'



  end subroutine transformBC

end module wedge_class
