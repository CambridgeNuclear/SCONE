module squareCylinder_class

  use numPrecision
  use universalVariables
  use genericProcedures,  only : fatalError, numToChar, swap
  use dictionary_class,   only : dictionary
  use surface_inter,      only : surface, kill_super => kill

  implicit none
  private

  !!
  !! Square cylinder aligned with x,y or z axis
  !!
  !! F(r) = maxval(abs(r - o) - a)
  !!
  !! Where: a -> 2D halfwidth vector, o-> 2D origin position
  !!        maxval -> maximum element (L_inf norm)
  !!
  !! Three diffrent types are avaliable
  !!   xSquareCylinder -> aligned with X-axis
  !!   ySquareCylinder -> aligned with Y-axis
  !!   zSquareCylinder -> aligned with Z-axis
  !!
  !! Surface Tolerance: SURF_TOL
  !!
  !! Sample Dictionary Input:
  !!   x { type xSquareCylinder; id 92; origin (0.0 0.0 9.0); halfwidth(0.0 2.0 0.3);}
  !!   y { type ySquareCylinder; id 92; origin (0.0 0.0 9.0); halfwidth(2.0 0.0 0.3);}
  !!
  !!   Halfwidth and origin entry in axis along the cylinder direction is ignored.
  !!
  !! Boundary Conditions:
  !!   BC order: x_min, x_max, y_min, y_max, z_min, z_max
  !!
  !!   Each face can have diffrent BC. Any combination is supported with co-ordinate transform.
  !!   BCs on all faces (even the infinate ones) must be provided. Of course BC for planes normal
  !!   to the cylinder axis does not matter.
  !!
  !! Private Members:
  !!   origin    -> position of the middle of the squareCylinder
  !!   halfwidth -> Halfwidths (half-length) of the squareCylinder in each direction (must be > 0.0)
  !!   BC        -> Boundary conditions flags for each face (x_min, x_max, y_min, y_max, z_min, z_max)
  !!   plane     -> Indexes of the axis, which are the plane of the cylinder
  !!       (e.g. for xSquareCylinder, y-z)
  !!   axis      -> Index of the axis of the plane
  !!
  !! Interface:
  !!   surface interface
  !!
  type, public, extends(surface) :: squareCylinder
    private
    real(defReal), dimension(2)     :: origin    = ZERO
    real(defReal), dimension(2)     :: halfwidth = ZERO
    integer(shortInt), dimension(6) :: BC = VACUUM_BC
    integer(shortInt), dimension(2) :: plane = 0
    integer(shortInt)               :: axis  = 0

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
    procedure :: explicitRayBC
    procedure :: transformBC
  end type squareCylinder


contains

  !!
  !! Return surface type name
  !!
  !! See surface_inter for more details
  !!
  pure function myType(self) result(str)
    class(squareCylinder), intent(in) :: self
    character(:), allocatable         :: str
    character(100), parameter :: Here = 'myType (squareCylinder_class.f90)'

    select case (self % axis)
      case (X_AXIS)
        str = 'xSquareCylinder'

      case (Y_AXIS)
        str = 'ySquareCylinder'

      case (Z_AXIS)
        str = 'zSquareCylinder'

      case default
        str = 'Unknown SquareCylinder'

    end select

  end function myType

  !!
  !! Initialise squareCylinder from a dictionary
  !!
  !! See surface_inter for more details
  !!
  !! Errors:
  !!  fatalError for -ve id or meningful halfspace. Non-meaningful halfspace entry is along
  !!    the culinder axis.
  !!
  subroutine init(self, dict)
    class(squareCylinder), intent(inout)     :: self
    class(dictionary), intent(in)            :: dict
    integer(shortInt)                        :: id, N
    real(defReal), dimension(:), allocatable :: temp
    character(nameLen)                       :: type
    character(100), parameter :: Here = 'init (squareCylinder_class.f90)'

    ! Load id
    call dict % get(id,'id')
    if (id <= 0) call fatalError(Here, 'ID must be <= 0. Is: '//numToChar(id))
    call self % setID(id)

    ! Select type
    call dict % get(type,'type')
    select case(type)
      case('xSquareCylinder')
        self % axis = X_AXIS
        self % plane = [Y_AXIS, Z_AXIS]

      case('ySquareCylinder')
        self % axis = Y_AXIS
        self % plane = [X_AXIS, Z_AXIS]

      case('zSquareCylinder')
        self % axis = Z_AXIS
        self % plane = [X_AXIS, Y_AXIS]

      case default
        call fatalError(Here, 'Unknown type of square cylinder: '//type)

    end select

    ! Load origin
    call dict % get(temp,'origin')
    N = size(temp)
    if (N /= 3) call fatalError(Here,'origin must have size 3. Has: '//numToChar(N))
    self % origin = temp(self % plane)


    ! Load halfwidth
    call dict % get(temp,'halfwidth')
    N = size(temp)
    if (N /= 3) then
      call fatalError(Here, 'halfwidth must have size 3. Has: '//numToChar(N))

    else if (any(temp(self % plane) < ZERO)) then
      call fatalError(Here, 'halfwidth cannot have -ve values.')

    end if
    self % halfwidth = temp(self % plane)

  end subroutine init

  !!
  !! Return axis-aligned bounding squareCylinder for the surface
  !!
  !! See surface_inter for details
  !!
  pure function boundingBox(self) result(aabb)
    class(squareCylinder), intent(in) :: self
    real(defReal), dimension(6)       :: aabb

    aabb(self % plane)     = self % origin - self % halfwidth
    aabb(self % plane + 3) = self % origin + self % halfwidth
    aabb(self % axis)      = -INF
    aabb(self % axis +3)   = INF

  end function boundingBox

  !!
  !! Evaluate surface expression c = F(r)
  !!
  !! See surface_inter for details
  !!
  pure function evaluate(self, r) result(c)
    class(squareCylinder), intent(in)       :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: c
    real(defReal), dimension(2)             :: rb

    ! Move to origin-frame and evaluate
    rb = r(self % plane) - self % origin
    c = maxval(abs(rb) - self % halfwidth)

  end function evaluate

  !!
  !! Return distance to the surface
  !!
  !! See surface_inter for details
  !!
  !! Uses algorithim similar to box_class.
  !!
  !! NOTE: Currently if a ray is PERFECTLY (to the bit) in one of the boundary planes
  !!   its NOT contained between planes ( allowed d is empty set). See note in box_class.
  !!
  pure function distance(self, r, u) result(d)
    class(squareCylinder), intent(in)       :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    real(defReal)                           :: d
    real(defReal), dimension(2)             :: rl, ul, rb, a_near, a_far
    real(defReal)                           :: far, near, test_far, test_near
    integer(shortInt)                       :: ax, i
    real(defReal), parameter :: FP_MISS_TOL = ONE + 10.0_defReal * epsilon(ONE)

    ! Get position and direction in the plane
    rl = r(self % plane)
    ul = u(self % plane)

    ! Transfrom to frame centered at the slab and choose
    ! nearest and furthest squareCylinder vertex
    rb = rl - self % origin
    a_far = sign(self % halfwidth, ul)
    a_near = -a_far

    ! Initialise intersection set
    far  = huge(far)
    near = -huge(near)

    ! Obtain set of distances inside the squareCylinder via set intersection
    ! Note that a_near/far; rb is 2D and u is 3D hence different indexing
    do i = 1, 2
      ax = self % plane(i)
      if (u(ax) /= ZERO ) then
        test_near = (a_near(i) - rb(i)) / u(ax)
        test_far = (a_far(i) - rb(i)) / u(ax)

      else ! Line is either fully inside or fully outside
        test_near = sign(INF, (a_near(i) - rb(i)))
        test_far  = sign(INF, (a_far(i) - rb(i)))

        ! If direction is -0.0 (intead of just 0.0) Order may be
        ! wrong and it is necesary to swap
        if (test_near > test_far) call swap(test_near, test_far)

      end if
      far = min(far, test_far)
      near = max(near, test_near)

    end do

    ! Choose correct distance for different cases
    if (far <= near * FP_MISS_TOL) then ! There is no intersection
      d = INF

    else if ( abs(self % evaluate(r)) < self % surfTol()) then ! Point at the surface
      ! Choose distance with largest absolute value
      if (abs(far) >= abs(near)) then
        d = far
      else
        d = near
      end if

    else ! Normal hit. Closest distance
      if (near <= ZERO) then
        d = far
      else
        d = near
      end if

    end if

    ! Cap the distance
    if (d <= ZERO .or. d > INF) then
      d = INF
    end if

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
    class(squareCylinder), intent(in)       :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    logical(defBool)                        :: halfspace
    real(defReal), dimension(2)             :: rl, ul
    integer(shortInt)                       :: maxCom
    real(defReal)                           :: proj

    ! Get position in the plane & direction
    rl = r(self % plane) - self % origin
    ul = u(self % plane)

    ! Find index of maximum component to identify axis
    ! of surface normal
    maxCom = maxloc(abs(rl) - self % halfwidth, 1)

    ! Projection of direction on the normal
    proj = ul(maxCom) * sign(ONE, rl(maxCom))

    halfspace = proj > ZERO

    ! Parallel direction
    ! Need to use position to determine halfspace
    if (proj == ZERO) then
      halfspace = self % evaluate(r) >= ZERO
    end if

  end function going

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(squareCylinder), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Local
    self % origin    = ZERO
    self % halfwidth = ZERO
    self % BC        = VACUUM_BC
    self % plane     = 0

  end subroutine kill

  !!
  !! Set boundary conditions
  !!
  !! See surface_inter for details
  !!
  subroutine setBC(self, BC)
    class(squareCylinder), intent(inout)        :: self
    integer(shortInt), dimension(:), intent(in) :: BC
    integer(shortInt)                           :: i
    character(100),parameter :: Here = 'setBC (squareCylinder_inter.f90)'

    if(size(BC) < 6) call fatalError(Here,'Wrong size of BC string. Must be at least 6')

    ! Load BC codes
    self % BC = BC(1:6)

    ! Verify that all BC flags make sense
    do i = 1, 6
      select case(BC(i))
        case (VACUUM_BC, REFLECTIVE_BC, PERIODIC_BC)
          ! Do nothing, pass
        case default
          call fatalError(Here,'Unrecognised BC: '//numToChar(BC(i))//' in position: '//numToChar(i))

      end select
    end do

    ! Verify periodic BCs
    if(.not.all( (self % BC([1,3,5] ) == PERIODIC_BC) .eqv. (self % BC([2,4,6]) == PERIODIC_BC))) then
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
    class(squareCylinder), intent(in)          :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u
    integer(shortInt)                          :: ax, bc, i
    real(defReal)                              :: r0
    character(100), parameter :: Here = 'explicitBC (squareCylinder_class.f90)'

    ! Loop over directions
    ! Becouse of the mix of 2D and 3D vectors to get right component use:
    !   i -> for 2D vectors
    !   ax -> for 3D vectors (r & u)
    axis : do i = 1, 2
      ax = self % plane(i)

      ! Find position wrt origin
      r0 = r(ax) - self % origin(i)

      ! Skip if particle is well inside the domain
      if (abs(r0)  <= self % halfwidth(i) - self % surfTol()) cycle axis

      ! Choose correct BC
      if (r0 < ZERO) then
        bc = self % BC(2*ax - 1)
      else
        bc = self % BC(2*ax)
      end if

      ! Apply BC
      select case(bc)
        case (VACUUM_BC)
          ! Do nothing. Pass

        case (REFLECTIVE_BC)
          u(ax) = -u(ax)

        case (PERIODIC_BC)
          ! Calculate displacement and perform translation
          r(ax) = r(ax) - TWO * sign(self % halfwidth(i), r0)

        case default
          call fatalError(Here, 'Unrecognised BC: '// numToChar(bc))

      end select
    end do axis

  end subroutine explicitBC

  !!
  !! Apply explicit BCs for ray problems: enforces a reflection when
  !! a vacuum was hit and reports this
  !!
  !! See surface_inter for details
  !!
  !! Note:
  !!   - Go through all directions in order to account for corners
  !!
  subroutine explicitRayBC(self, r, u, hitVacuum)
    class(squareCylinder), intent(in)          :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u
    logical(defBool), intent(out)              :: hitVacuum
    integer(shortInt)                          :: ax, bc, i
    real(defReal)                              :: r0
    character(100), parameter :: Here = 'explicitRayBC (squareCylinder_class.f90)'

    hitVacuum = .FALSE.

    ! Loop over directions
    ! Becouse of the mix of 2D and 3D vectors to get right component use:
    !   i -> for 2D vectors
    !   ax -> for 3D vectors (r & u)
    axis : do i = 1, 2
      ax = self % plane(i)

      ! Find position wrt origin
      r0 = r(ax) - self % origin(i)

      ! Skip if particle is well inside the domain
      if (abs(r0)  <= self % halfwidth(i) - self % surfTol()) cycle axis

      ! Choose correct BC
      if (r0 < ZERO) then
        bc = self % BC(2*ax - 1)
      else
        bc = self % BC(2*ax)
      end if

      ! Apply BC
      select case(bc)
        case (VACUUM_BC)
          ! Treat as reflective but state that a vacuum was struck
          u(ax) = -u(ax)
          hitVacuum = .TRUE.

        case (REFLECTIVE_BC)
          u(ax) = -u(ax)

        case (PERIODIC_BC)
          ! Calculate displacement and perform translation
          r(ax) = r(ax) - TWO * sign(self % halfwidth(i), r0)

        case default
          call fatalError(Here, 'Unrecognised BC: '// numToChar(bc))

      end select
    end do axis

  end subroutine explicitRayBC

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
    class(squareCylinder), intent(in)          :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u
    real(defReal), dimension(2)                :: a_bar
    integer(shortInt), dimension(2)            :: Ri
    integer(shortInt)                          :: ax, t, bc, i
    real(defReal)                              :: a0, d, r0
    character(100), parameter :: Here = 'transformBC (squareCylinder_class.f90)'

    ! Calculate halfwidth reduced by the surface_tolerance
    ! Necessary to capture particles at the boundary
    a_bar = self % halfwidth - self % surfTol()

    ! Calculate distance (in # of transformations) in each direction
    Ri = ceiling(abs(r(self % plane) - self % origin) / a_bar) / 2

    ! Loop over directions & number of transformations
    ! Becouse of the mix of 2D and 3D vectors to get right component use:
    !   i -> for 2D vectors
    !   ax -> for 3D vectors (r & u)
    axis : do i = 1, 2
      ax = self % plane(i)

      trans : do t = 1, Ri(i)
        ! Find position wrt origin
        r0 = r(ax) - self % origin(i)

        ! Choose correct BC
        if ( r0 < ZERO) then
          bc = self % BC(2*ax - 1)
        else
          bc = self % BC(2*ax)
        end if

        ! Apply BC
        select case(bc)
          case (VACUUM_BC)
            ! Do nothing. Pass

          case (REFLECTIVE_BC)
            ! Find position of the plane
            a0 = sign(self % halfwidth(i), r0) + self % origin(i)

            ! Calculate displacment and perform reflection
            d = r(ax) - a0
            r(ax) = r(ax) - TWO * d
            u(ax) = -u(ax)

          case (PERIODIC_BC)
            ! Calculate displacement and perform translation
            d = sign(self % halfwidth(i), r0)
            r(ax) = r(ax) - TWO * d

          case default
            call fatalError(Here, 'Unrecognised BC: '// numToChar(bc))
        end select

      end do trans
    end do axis

  end subroutine transformBC

end module squareCylinder_class
