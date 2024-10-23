module truncCylinder_class

  use numPrecision
  use universalVariables
  use genericProcedures,  only : fatalError, numToChar, swap, isEqual
  use dictionary_class,   only : dictionary
  use surface_inter,      only : surface, kill_super => kill

  implicit none
  private

  !!
  !! Finite length cylinder aligned with one of the co-ord axis (x,y or z)
  !!
  !! F(r) = max[ {(r1-o1)^2 + (r2-o2)^2 - R^2}/(2R); abs(r3-o3)-a]
  !!
  !! r -> position; o -> origin; R -> radius; a-> halfwidth
  !! 1,2 -> planar axis, 3-> cylinder axis
  !! Denominator 2*R is required for consistant surface tolerance thickness
  !!
  !! Three diffrent types are avaliable
  !!   xTruncCylinder -> aligned with X-axis
  !!   yTruncCylinder -> aligned with Y-axis
  !!   zTruncCylinder -> aligned with Z-axis
  !!
  !! Surface Tolerance: SURF_TOL
  !!
  !! Sample Dictionary Input:
  !!  x {type xTruncCylinder; id 2; origin (1.0 -2.0 0.0); halfwidth 1.2; radius 0.5;}
  !!  y {type yTruncCylinder; id 2; origin (1.0 2.0 7.0); halfwidth 1.3; radius 1.5;}
  !!
  !! Boundary Conditions:
  !!   BC order: a_min, a_max
  !!   Where a is cylinder axis [x,y,z]
  !!   BC in radial direction of the cylinder is always VACUUM
  !!
  !! Private Members:
  !!   origin -> Position of the centre of the cylinder
  !!   a      -> Axial halfwidth (>0.0)
  !!   r      -> Radius (>0.0)
  !!   axis   -> Cylinder axis specifier
  !!   plane  -> Planar axis specifiers
  !!   BC     -> Boundary conditions flags [a_min, a_max]
  !!
  !! Interface:
  !!   surface interface
  !!
  type, public, extends(surface) :: truncCylinder
    private
    real(defReal), dimension(3)     :: origin = ZERO
    real(defReal)                   :: a      = ZERO
    real(defReal)                   :: r      = ZERO
    integer(shortInt)               :: axis   = -7
    integer(shortInt), dimension(2) :: plane  = -7
    integer(shortInt), dimension(2) :: BC     = VACUUM_BC

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

  end type truncCylinder

contains

  !!
  !! Return surface type name
  !!
  !! See surface_inter for more details
  !!
  pure function myType(self) result(str)
    class(truncCylinder), intent(in) :: self
    character(:), allocatable         :: str
    character(100), parameter :: Here = 'myType (truncCylinder_class.f90)'

    select case (self % axis)
      case (X_AXIS)
        str = 'xTruncCylinder'

      case (Y_AXIS)
        str = 'yTruncCylinder'

      case (Z_AXIS)
        str = 'zTruncCylinder'

      case default
        str = 'Unknown Truncated Cylinder'

    end select

  end function myType

  !!
  !! Initialise truncCylinder from a dictionary
  !!
  !! See surface_inter for more details
  !!
  subroutine init(self, dict)
    class(truncCylinder), intent(inout)      :: self
    class(dictionary), intent(in)            :: dict
    integer(shortInt)                        :: id
    real(defReal), dimension(:), allocatable :: temp
    real(defReal)                            :: r, a
    character(nameLen)                       :: type
    character(100), parameter :: Here = 'init (truncCylinder_class.f90)'

    ! Read id
    call dict % get(id, 'id')
    if (id <= 0) call fatalError(Here,'ID must be <= 0. Is: '//numToChar(id))
    call self % setID(id)

    ! Read origin
    call dict % get(temp, 'origin')
    if (size(temp) /= 3) then
      call fatalError(Here, 'origin must have size 3. Has: '//numToChar(size(temp)))
    end if
    self % origin = temp

    ! Read radius
    call dict % get(r, 'radius')
    if (r <= ZERO) call fatalError(Here, 'Radius must be +ve. Is: '//numToChar(r))
    self % r = r

    ! Read halfwidth
    call dict % get(a, 'halfwidth')
    if (a <= ZERO) call fatalError(Here, 'Halfwidth must be +ve. Is: '//numToChar(a))
    self % a = a

    ! Read type
    call dict % get(type, 'type')
    select case (type)
      case ('xTruncCylinder')
        self % axis = X_AXIS
        self % plane = [Y_AXIS, Z_AXIS]

      case ('yTruncCylinder')
        self % axis = Y_AXIS
        self % plane = [X_AXIS, Z_AXIS]

      case ('zTruncCylinder')
        self % axis = Z_AXIS
        self % plane = [X_AXIS, Y_AXIS]

      case default
        call fatalError(Here, 'Unknown type of truncCylindeer: '//type)

    end select

  end subroutine init

  !!
  !! Return axis-aligned bounding truncCylinder for the surface
  !!
  !! See surface_inter for details
  !!
  pure function boundingBox(self) result(aabb)
    class(truncCylinder), intent(in) :: self
    real(defReal), dimension(6)      :: aabb

    ! Set values along the axis
    aabb(self % axis) = self % origin(self % axis) - self % a
    aabb(self % axis + 3) = self % origin(self % axis) + self % a

    ! Set values in plane
    aabb(self % plane) = self % origin(self % plane) - self % r
    aabb(self % plane + 3) = self % origin(self % plane) + self % r

  end function boundingBox

  !!
  !! Evaluate surface expression c = F(r)
  !!
  !! See surface_inter for details
  !!
  pure function evaluate(self, r) result(c)
    class(truncCylinder), intent(in)       :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: c
    integer(shortInt), dimension(2)         :: p
    integer(shortInt)                       :: a

    ! Get plane and axis indexes into shorter variables (for clarity)
    p = self % plane
    a = self % axis

    ! Evaluate surface expression
    c = (sum((r(p) - self % origin(p))**2) - self % r*self % r) / self % r * HALF
    c = max(c, abs(r(a) - self % origin(a)) - self % a)

  end function evaluate

  !!
  !! Return distance to the surface
  !!
  !! See surface_inter for details
  !!
  !! Uses intersection method analogous to box_class
  !! For min and max distance calculation for cylinder see cylinder_class
  !!
  pure function distance(self, r, u) result(d)
    class(truncCylinder), intent(in)        :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    real(defReal)                           :: d
    integer(shortInt), dimension(2)         :: p
    integer(shortInt)                       :: ax
    real(defReal)                           :: c1, k, delta, a, near, far
    real(defReal)                           :: test_near, test_far, rb, c
    real(defReal), parameter :: FP_MISS_TOL = ONE + 10.0_defReal * epsilon(ONE)

    ! Get plane and axis indexes into shorter variables (for clarity)
    p = self % plane
    ax = self % axis

    ! Parameters of quadratic expression
    c1 = sum((r(p) - self % origin(p))**2) - self % r * self % r
    k = dot_product(r(p) - self % origin(p) , u(p))
    a = ONE - u(ax)**2
    delta = k*k - a*c1  ! Technically delta/4

    ! Find closes & furthest distance
    if (delta <= ZERO .or. isEqual(a, ZERO)) then ! No intersection
      far = INF
      near = sign(INF, c1) ! If ray is parallel inside the cylinder it must be fully contained

    else
      far = (-k + sqrt(delta)) / a
      near = (-k - sqrt(delta)) / a
      if (far < near) call swap(far, near) ! Ensure correct order for any orientation

    end if

    ! Axial distance
    rb = r(ax) - self % origin(ax)
    if (u(ax) /= ZERO) then
      test_near = (-self % a - rb) / u(ax)
      test_far = (self % a - rb) / u(ax)

    else
      test_near = sign(INF, (-self % a - rb))
      test_far = sign(INF, (self % a - rb))

    end if

    ! Ensure correct order for any orientation
    if (test_far < test_near) call swap(test_far, test_near)

    ! Get intersection
    far = min(far, test_far)
    near = max(near, test_near)

    ! Evaluate surface expression
    c = max(c1 / self % r * HALF, abs(rb) - self % a)

    ! Choose correct distance for different cases
    if (far <= near * FP_MISS_TOL) then ! There is no intersection
      d = INF

    else if ( abs(c) < self % surfTol()) then ! Point at the surface
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
  !! Selects closest face (by normal distance to the face -> value of c in surf-expresion).
  !! Then it uses the projection of direction on the normal of that face to select next
  !! halfspace.
  !!
  pure function going(self, r, u) result(halfspace)
    class(truncCylinder), intent(in)       :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    logical(defBool)                        :: halfspace
    real(defReal)                           :: c1, c2 , rv, proj, c
    real(defReal), dimension(2)             :: rp
    integer(shortInt), dimension(2)         :: p
    integer(shortInt)                       :: a

    ! Get plane and axis indexes into shorter variables (for clarity)
    p = self % plane
    a = self % axis

    ! Determine closest surface
    rp = r(p) - self % origin(p)
    c1 = (sum(rp**2) - self % r * self % r) / self % r * HALF

    rv = r(a) - self % origin(a)
    c2 = abs(rv) - self % a

    if (c1 >= TWO * self % r * c2) then
      proj = dot_product(u(p), rp)
      c = c1

    else
      proj = u(a) * rv
      c = c2

    end if

    ! Parallel direction. Need to use position to determine halfspace.
    if (isEqual(proj, ZERO)) then
      halfspace = c >= ZERO
      return
      
    end if

    ! Determine next halfspace
    halfspace = proj > ZERO

  end function going

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(truncCylinder), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Local
    self % origin = ZERO
    self % a      = ZERO
    self % r      = ZERO
    self % axis   = -7
    self % plane  = -7
    self % BC     = VACUUM_BC

  end subroutine kill

  !!
  !! Set boundary conditions
  !!
  !! See surface_inter for details
  !!
  subroutine setBC(self, BC)
    class(truncCylinder), intent(inout)         :: self
    integer(shortInt), dimension(:), intent(in) :: BC
    integer(shortInt)                           :: i
    character(100),parameter :: Here = 'setBC (truncCylinder_inter.f90)'

    if(size(BC) < 2) call fatalError(Here,'Wrong size of BC string. Must be at least 2')

    ! Load BC codes
    self % BC = BC(1:2)

    ! Verify that all BC flags make sense
    do i = 1, 2
      select case(BC(i))
        case (VACUUM_BC, REFLECTIVE_BC, PERIODIC_BC)
          ! Do nothing, pass
        case default
          call fatalError(Here,'Unrecognised BC: '//numToChar(BC(i))//' in position: '//numToChar(i))

      end select
    end do

    ! Verify periodic BCs
    if((self % BC(1) == PERIODIC_BC) .neqv. (self % BC(1) == PERIODIC_BC)) then
      call fatalError(Here,'Periodic BC need to be applied to opposite surfaces')

    end if

  end subroutine setBC

  !!
  !! Apply explicit BCs
  !!
  !! See surface_inter for details
  !!
  !! Approach analogous to box_class applied to axial direction only
  !!
  subroutine explicitBC(self, r, u)
    class(truncCylinder), intent(in)           :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u
    integer(shortInt)                          :: ax, bc
    real(defReal)                              :: r0
    character(100), parameter :: Here = 'explicitBC (truncCylinder_class.f90)'

    ax = self % axis

    ! Find position wrt origin
    r0 = r(ax) - self % origin(ax)

    ! Do nothing if particle is well inside the domain
    if (abs(r0)  <= self % a - self % surfTol()) return

    ! Choose correct BC
    if (r0 < ZERO) then
      bc = self % BC(1)
    else
      bc = self % BC(2)
    end if

    ! Apply BC
    select case(bc)
      case (VACUUM_BC)
        ! Do nothing. Pass

      case (REFLECTIVE_BC)
        u(ax) = -u(ax)

      case (PERIODIC_BC)
        ! Calculate displacement and perform translation
        r(ax) = r(ax) - TWO * sign(self % a, r0)

      case default
        call fatalError(Here, 'Unrecognised BC: '// numToChar(bc))

    end select

  end subroutine explicitBC

  !!
  !! Apply co-ordinate transform BC
  !!
  !! See surface_inter for details
  !!
  !! Approach analogous to box_class applied to axial direction only
  !!
  subroutine transformBC(self, r, u)
    class(truncCylinder), intent(in)           :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u
    integer(shortInt)                          :: Ri, ax, t, bc
    real(defReal)                              :: a_bar, a0, d, r0
    character(100), parameter :: Here = 'transformBC (truncCylinder_class.f90)'

    ! Select axis
    ax = self % axis

    ! Calculate halfwidth reduced by the surface_tolerance
    ! Necessary to capture particles at the boundary
    a_bar = self % a - self % surfTol()

    ! Calculate distance (in # of transformations) in each direction
    Ri = ceiling(abs(r(ax) - self % origin(ax)) / a_bar) / 2

    trans : do t = 1, Ri
      ! Find position wrt origin
      r0 = r(ax) - self % origin(ax)

      ! Choose correct BC
      if ( r0 < ZERO) then
        bc = self % BC(1)
      else
        bc = self % BC(2)
      end if

      ! Apply BC
      select case(bc)
        case (VACUUM_BC)
          ! Do nothing. Pass

        case (REFLECTIVE_BC)
          ! Find position of the plane
          a0 = sign(self % a, r0) + self % origin(ax)

          ! Calculate displacment and perform reflection
          d = r(ax) - a0
          r(ax) = r(ax) - TWO * d
          u(ax) = -u(ax)

        case (PERIODIC_BC)
          ! Calculate displacement and perform translation
          d = sign(self % a, r0)
          r(ax) = r(ax) - TWO * d

        case default
          call fatalError(Here, 'Unrecognised BC: '// numToChar(bc))
      end select
    end do trans

  end subroutine transformBC

end module truncCylinder_class
