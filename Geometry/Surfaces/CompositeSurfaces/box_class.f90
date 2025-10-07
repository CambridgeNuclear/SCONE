module box_class

  use numPrecision
  use universalVariables
  use genericProcedures,  only : fatalError, numToChar, swap
  use dictionary_class,   only : dictionary
  use compSurface_inter,  only : compSurface
  use surface_inter,      only : kill_super => kill

  implicit none
  private

  !!
  !! Axis Aligned Box
  !!
  !! F(r) = maxval(abs(r - o) - a)
  !!
  !! Where: a -> halfwidth vector, o -> origin position
  !!        maxval -> maximum element (L_inf norm)
  !!
  !! Surface Tolerance: SURF_TOL
  !!
  !! Sample Dictionary Input:
  !!   aab { type box; id 92; origin (0.0 0.0 9.0); halfwidth(1.0 2.0 0.3);}
  !!
  !! Boundary Conditions:
  !!   BC order: x_min, x_max, y_min, y_max, z_min, z_max
  !!   Each face can have diffrent BC. Any combination is supported with co-ordinate transform.
  !!
  !! Private Members:
  !!   origin -> poosition of the middle of the box
  !!   halfwidth -> Halfwidths (half-length) of the box in each direction (must be > 0.0)
  !!   BC -> Boundary conditions flags for each face (x_min, x_max, y_min, y_max, z_min, z_max)
  !!
  !! Interface:
  !!   surface interface
  !!
  type, public, extends(compSurface) :: box
    private
    real(defReal), dimension(3)     :: origin    = ZERO
    real(defReal), dimension(3)     :: halfwidth = ZERO
    integer(shortInt), dimension(6) :: BC = VACUUM_BC

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
  end type box


contains

  !!
  !! Return surface type name
  !!
  !! See surface_inter for more details
  !!
  pure function myType(self) result(str)
    class(box), intent(in) :: self
    character(:), allocatable  :: str

    str = 'box'

  end function myType

  !!
  !! Initialise box from a dictionary
  !!
  !! See surface_inter for more details
  !!
  !! Errors:
  !!
  subroutine init(self, dict)
    class(box), intent(inout)     :: self
    class(dictionary), intent(in) :: dict
    integer(shortInt)                        :: id, N
    real(defReal), dimension(:), allocatable :: temp
    character(100), parameter :: Here = 'init (box_class.f90)'

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
    call dict % get(temp,'halfwidth')
    N = size(temp)
    if (N /= 3) then
      call fatalError(Here, 'halfwidth must have size 3. Has: '//numToChar(N))

    else if (any(temp < ZERO)) then
      call fatalError(Here, 'halfwidth cannot have -ve values.')

    end if
    self % halfwidth = temp

  end subroutine init

  !!
  !! Return axis-aligned bounding box for the surface
  !!
  !! See surface_inter for details
  !!
  pure function boundingBox(self) result(aabb)
    class(box), intent(in)      :: self
    real(defReal), dimension(6) :: aabb

    aabb(1:3) = self % origin - self % halfwidth
    aabb(4:6) = self % origin + self % halfwidth

  end function boundingBox

  !!
  !! Evaluate surface expression c = F(r)
  !!
  !! See surface_inter for details
  !!
  pure function evaluate(self, r) result(c)
    class(box), intent(in)                  :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: c
    real(defReal), dimension(3)             :: rb

    ! Move to origin-frame and evaluate
    rb = r - self % origin
    c = maxval(abs(rb) - self % halfwidth)

  end function evaluate

  !!
  !! Return distance to the surface
  !!
  !! See surface_inter for details
  !!
  !! Distance calculation works by:
  !!   1) Note that the box is intersection of space between parallel planes
  !!   2) For each pair of planes calculate the section of the distance line, which lies
  !!      between the planes
  !!   3) Find intersection of these sections (by keeping maximum near distance and minimum
  !!      far distance)
  !!   4) If set is empty (far <= near) {well... has measure 0.0 technically} ray is a miss.
  !!      Otherwise either near or far will be the closest distance depending on whether
  !!      the particle is in or outside the box
  !!
  !! NOTE: The algorithim as implemented may have some problems in the corners. If particle
  !!  is EXACTLY (to the bit) in the corner, or it goes through corner while being PERFECTLY
  !!  (again to the bit) in one of the planes
  !!
  !!  The reason for this is that a ray contained EXACTLY in one of the boundary planes
  !!  may be assigned with an empty section (instead of infinate length).
  !!
  pure function distance(self, r, u) result(d)
    class(box), intent(in)              :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    real(defReal)                           :: d
    real(defReal), dimension(3)             :: rb, a_near, a_far
    real(defReal)                           :: far, near, test_far, test_near
    integer(shortInt)                       :: ax
    real(defReal), parameter :: FP_MISS_TOL = ONE + 10.0_defReal * epsilon(ONE)

    ! Transfrom to frame centered at the slab and choose
    ! nearest and furthest box vertex
    rb = r - self % origin
    a_far = sign(self % halfwidth, u)
    a_near = -a_far

    ! Initialise intersection set
    far  = huge(far)
    near = -huge(near)

    ! Obtain set of distances inside the box via set intersection
    do ax = 1, 3
      if (u(ax) /= ZERO ) then
        test_near = (a_near(ax) - rb(ax)) / u(ax)
        test_far = (a_far(ax) - rb(ax)) / u(ax)

      else  ! Line is either fully inside or fully outside

        ! TODO: In current setup line may not inclusive of boundary planes i.e.
        !   if either (a_near - r) or (a_far - r) is 0.0 then both INF may have the same sign
        test_near = sign(INF, (a_near(ax) - rb(ax)))
        test_far  = sign(INF, (a_far(ax) - rb(ax)))

        ! If direction is -0.0 (intead of just 0.0) Order may be
        ! wrong and it is necesaryto swap
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
    class(box), intent(in)              :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    logical(defBool)                        :: halfspace
    real(defReal), dimension(3)             :: rl
    integer(shortInt)                       :: maxCom
    real(defReal)                           :: proj

    rl = r - self % origin

    ! Find index of maximum component to identify axis
    ! of surface normal
    maxCom = maxloc(abs(rl) - self % halfwidth, 1)

    ! Projection of direction on the normal
    proj = u(maxCom) * sign(ONE, rl(maxCom))

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
    class(box), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Local
    self % origin = ZERO
    self % halfwidth = ZERO
    self % BC = VACUUM_BC

  end subroutine kill

  !!
  !! Set boundary conditions
  !!
  !! See surface_inter for details
  !!
  subroutine setBC(self, BC)
    class(box), intent(inout)                   :: self
    integer(shortInt), dimension(:), intent(in) :: BC
    integer(shortInt)                           :: i
    character(100),parameter :: Here = 'setBC (box_class.f90)'

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
    if (.not. all( (self % BC([1,3,5]) == PERIODIC_BC) .eqv. (self % BC([2,4,6]) == PERIODIC_BC))) then
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
    class(box), intent(in)                     :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u
    integer(shortInt)                          :: ax, bc
    real(defReal)                              :: r0
    character(100), parameter :: Here = 'explicitBC (box_class.f90)'

    ! Loop over directions
    axis : do ax = 1, 3
      ! Find position wrt origin
      r0 = r(ax) - self % origin(ax)

      ! Skip if particle is well inside the domain
      if (abs(r0)  <= self % halfwidth(ax) - self % surfTol()) cycle axis

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
          r(ax) = r(ax) - TWO * sign(self % halfwidth(ax), r0)

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
    class(box), intent(in)                     :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u
    real(defReal), dimension(3)                :: a_bar
    integer(shortInt), dimension(3)            :: Ri
    integer(shortInt)                          :: ax, t, bc
    real(defReal)                              :: a0, d, r0
    character(100), parameter :: Here = 'transformBC (box_class.f90)'

    ! Calculate halfwidth reduced by the surface_tolerance
    ! Necessary to capture particles at the boundary
    a_bar = self % halfwidth - self % surfTol()

    ! Calculate distance (in # of transformations) in each direction
    Ri = ceiling(abs(r - self % origin) / a_bar) / 2

    ! Loop over directions & number of transformations
    axis : do ax = 1, 3
      trans : do t = 1, Ri(ax)
        ! Find position wrt origin
        r0 = r(ax) - self % origin(ax)

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
            a0 = sign(self % halfwidth(ax), r0) + self % origin(ax)

            ! Calculate displacment and perform reflection
            d = r(ax) - a0
            r(ax) = r(ax) - TWO * d
            u(ax) = -u(ax)

          case (PERIODIC_BC)
            ! Calculate displacement and perform translation
            d = sign(self % halfwidth(ax), r0)
            r(ax) = r(ax) - TWO * d

          case default
            call fatalError(Here, 'Unrecognised BC: '// numToChar(bc))
        end select

      end do trans
    end do axis

  end subroutine transformBC

end module box_class
