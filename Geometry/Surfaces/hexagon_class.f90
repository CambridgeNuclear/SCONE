module hexagon_class

  use numPrecision
  use universalVariables
  use genericProcedures,  only : fatalError, numToChar, swap
  use dictionary_class,   only : dictionary
  use surface_inter,      only : surface, kill_super => kill

  implicit none
  private

  integer(shortInt), parameter :: pointType = 1, flatType = 2

  !!
  !! Hexagon aligned with x, y or z axis.
  !! The hexagon can also be oriented either to have the point or flat
  !! aligning with the first of the remaining two axis.
  !!
  !! Three different axis types are avaliable
  !!   xHexagon -> aligned with X-axis
  !!   yHexagon -> aligned with Y-axis
  !!   zHexagon -> aligned with Z-axis
  !!
  !! One can then specify the orientation to be type 1 or 2.
  !! Type 1 has the flat surface perpendicular to the first of the two axis in the plane
  !! Type 2 has the flat surface parallel to the first of the two axis in the plane
  !!
  !! Surface Tolerance: SURF_TOL
  !!
  !! Sample Dictionary Input:
  !!   x { type xHexagon; id 92; orientation 1; origin (0.0 0.0 9.0); halfwidth 7;}
  !!   y { type yHexagon; id 94; orientation 2; origin (0.0 0.0 9.0); halfwidth 3;}
  !!
  !!  Origin entry in axis along the cylinder direction is ignored.
  !!  Halfwidth is the flat-to-flat distance.
  !!
  !! Boundary Conditions:
  !!   Only a single BC value.
  !!
  !!   For now, each face must have the same BC. Supports either periodic or vacuum BCs.
  !!
  !! Private Members:
  !!   origin      -> position of the middle of the hexagon
  !!   halfwidth   -> Halfwidth of the hexagon from one flat to the opposite. Same in all directions
  !!   BC          -> Boundary conditions flag for all faces in the plane
  !!   plane       -> Indices of the axis, which are the plane of the cylinder
  !!                  (e.g. for xHexagon, y-z)
  !!   axis        -> Index of the axis of the plane
  !!   flatAxis    -> Index into plane, deciding which plane axis has the flat of the hexagon
  !!   pointAxis   -> Index into plane, deciding which plane axis has the point of the hexagon
  !!
  !! Interface:
  !!   surface interface
  !!
  type, public, extends(surface) :: hexagon
    private
    real(defReal), dimension(2)     :: origin    = ZERO
    real(defReal)                   :: halfwidth = ZERO
    integer(shortInt)               :: BC = VACUUM_BC
    integer(shortInt), dimension(2) :: plane = 0
    integer(shortInt)               :: axis  = 0
    integer(shortInt)               :: flatAxis = 0
    integer(shortInt)               :: pointAxis = 0
    real(defReal), dimension(6,2)   :: verts = ZERO

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
  end type hexagon


contains

  !!
  !! Return surface type name
  !!
  !! See surface_inter for more details
  !!
  pure function myType(self) result(str)
    class(hexagon), intent(in) :: self
    character(:), allocatable  :: str
    character(100), parameter  :: Here = 'myType (hexagon_class.f90)'

    select case (self % axis)
      case (X_AXIS)
        str = 'xHexagon'

      case (Y_AXIS)
        str = 'yHexagon'

      case (Z_AXIS)
        str = 'zHexagon'

      case default
        str = 'Unknown hexagon'

    end select

  end function myType

  !!
  !! Initialise hexagon from a dictionary
  !!
  !! See surface_inter for more details
  !!
  !! Errors:
  !!  fatalError for -ve id or meaningful halfspace.
  !!  fatalError for unrecognised oritnetation
  !!
  subroutine init(self, dict)
    class(hexagon), intent(inout)            :: self
    class(dictionary), intent(in)            :: dict
    integer(shortInt)                        :: id, N, i, orient
    real(defReal), dimension(:), allocatable :: temp
    character(nameLen)                       :: type
    real(defReal), dimension(6)              :: angles
    character(100), parameter :: Here = 'init (hexagon_class.f90)'

    ! Load id
    call dict % get(id,'id')
    if (id <= 0) call fatalError(Here, 'ID must be <= 0. Is: '//numToChar(id))
    call self % setID(id)

    ! Select type
    call dict % get(type,'type')
    select case(type)
      case('xHexagon')
        self % axis = X_AXIS
        self % plane = [Y_AXIS, Z_AXIS]

      case('yHexagon')
        self % axis = Y_AXIS
        self % plane = [X_AXIS, Z_AXIS]

      case('zHexagon')
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
    call dict % get(self % halfwidth,'halfwidth')
    if (self % halfwidth < ZERO) call fatalError(Here, 'halfwidth cannot have a -ve value.')

    ! Get orientation
    call dict % get(orient, 'orientation')
    if (all(orient /= [pointType, flatType])) then
      call fatalError(Here, 'Unrecognised hexagon orientation: '//numToChar(orient))
    end if

    if (orient == pointType) then
      self % flatAxis = 2
      self % pointAxis = 1
      angles = [(i, i = 0,5)] * PI/3 + PI/6
    else
      self % flatAxis = 1
      self % pointAxis = 2
      angles = [(i, i = 0,5)] * PI/3
    end if
    self % verts(:,1) = TWO_SQRT3 * self % halfwidth * cos(angles) + self % origin(1)
    self % verts(:,2) = TWO_SQRT3 * self % halfwidth * sin(angles) + self % origin(2)

    
  end subroutine init

  !!
  !! Return axis-aligned bounding squareCylinder for the surface
  !!
  !! See surface_inter for details
  !!
  pure function boundingBox(self) result(aabb)
    class(hexagon), intent(in)  :: self
    real(defReal), dimension(6) :: aabb
    real(defReal)               :: originF, originP
    integer(shortInt)           :: f, p

    ! Two of the planes will be parallel to either the first or second axis in the plane
    f = self % plane(self % flatAxis)
    p = self % plane(self % pointAxis)
    originF = self % origin(self % flatAxis)
    originP = self % origin(self % pointAxis)

    aabb(f)     = originF - self % halfwidth
    aabb(f + 3) = originF + self % halfwidth
      
    ! The other two are shifted further by 2/sqrt(3)
    aabb(p)     = originP - self % halfwidth * TWO_SQRT3
    aabb(p + 3) = originP + self % halfwidth * TWO_SQRT3

    aabb(self % axis)      = -INF
    aabb(self % axis + 3)  = INF

  end function boundingBox

  !!
  !! Evaluate surface expression c = F(r)
  !!
  !! See surface_inter for details
  !!
  pure function evaluate(self, r) result(c)
    class(hexagon), intent(in)              :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: c
    real(defReal), dimension(2)             :: rb
    real(defReal)                           :: p, q, h

    ! Move to origin-frame and evaluate
    rb = abs(r(self % plane) - self % origin)
    
    p = rb(self % flatAxis)
    q = rb(self % pointAxis)

    h = self % halfwidth * TWO_SQRT3

    c = maxval([q - h * HALF, q + sqrt(3.0_defReal) * p - h])

  end function evaluate

  !!
  !! Return distance to the surface
  !!
  !! See surface_inter for details
  !!
  pure function distance(self, r, u) result(d)
    class(hexagon), intent(in)              :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    real(defReal)                           :: d
    real(defReal), dimension(2)             :: rl, ul, a1, a2, e, w
    real(defReal)                           :: d0, invDet
    integer(shortInt)                       :: i, iNext

    ! Get position and direction in the plane and centre
    rl = r(self % plane) - self % origin
    ul = u(self % plane)

    d = INF

    do i = 1, 6

      iNext = mod(i, 6) + 1
      
      ! Create the plane being intersected
      a1 = self % verts(i,:)
      a2 = self % verts(iNext,:)
      e = a2 - a1

      ! Ensure particle isn't running parallel
      invDet = ONE / (ul(1) * e(2) - ul(2) * e(1))
      if (abs(invDet) < 1E-10) cycle

      ! Invert the linear system
      w = rl - a1
      d0 = invDet * (e(2) * w(1) - e(1) * w(2))
      
      ! Ensure intersection is in front of ray
      if (d0 < ZERO .or. d0 > d) cycle

      d = d0

    end do

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
  !!   For parallel direction halfspace is assigned by the sign of `evaluate` result.
  !!
  pure function going(self, r, u) result(halfspace)
    class(hexagon), intent(in)              :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    logical(defBool)                        :: halfspace
    real(defReal), dimension(2)             :: rl, ul, a1, a2, e, n, w
    real(defReal)                           :: invDet, d, d0, proj
    integer(shortInt)                       :: i, iNext

    ! Get position in the plane & direction
    rl = r(self % plane) - self % origin
    ul = u(self % plane)

    d = INF
    do i = 1, 6
      
    iNext = mod(i, 6) + 1
      
      ! Create the plane being intersected
      a1 = self % verts(i,:)
      a2 = self % verts(iNext,:)
      e = a2 - a1

      ! Ensure particle isn't running parallel
      invDet = ONE / (ul(1) * e(2) - ul(2) * e(1))
      if (abs(invDet) < 1E-10) cycle

      ! Invert the linear system
      w = rl - a1
      d0 = invDet * (e(2) * w(1) - e(1) * w(2))
      
      ! Ensure intersection is in front of ray
      if (d0 < ZERO .or. d0 > d) cycle

      d = d0

      ! Projection of direction on the plane normal
      n = [-e(2), e(1)]
      n = n / norm2(n)
      proj = dot_product(ul, n)
      halfspace = proj > ZERO

    end do

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
    class(hexagon), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Local
    self % origin      = ZERO
    self % halfwidth   = ZERO
    self % BC          = VACUUM_BC
    self % plane       = 0
    self % flatAxis    = 0
    self % pointAxis   = 0
    self % verts       = ZERO

  end subroutine kill

  !!
  !! Set boundary conditions
  !!
  !! See surface_inter for details
  !!
  subroutine setBC(self, BC)
    class(hexagon), intent(inout)               :: self
    integer(shortInt), dimension(:), intent(in) :: BC
    integer(shortInt)                           :: i
    character(100),parameter                    :: Here = 'setBC (hexagon_class.f90)'

    if(size(BC) < 6) call fatalError(Here, 'Wrong size of BC string. Must be at least 6')
    if(all(BC /= BC(1))) call fatalError(Here, 'All BCs must be identical')

    ! Load BC codes
    self % BC = BC(1)

    ! Verify that all BC flags make sense
    select case(self % BC)
      case (VACUUM_BC, PERIODIC_BC)
          ! Do nothing, pass

      case (REFLECTIVE_BC)
        call fatalError(Here,'Reflective BCs not supported')

      case default
        call fatalError(Here,'Unrecognised BC: '//numToChar(BC))

    end select

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
    class(hexagon), intent(in)                 :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u
    integer(shortInt)                          :: bc, i, iHit, iNext
    real(defReal)                              :: r0, d, d0
    real(defReal), dimension(2)                :: rl, ul, e, b, proj, n
    character(100), parameter :: Here = 'explicitBC (hexagon_class.f90)'

    ! Apply BC
    select case(self % BC)
      case (VACUUM_BC)
        ! Do nothing. Pass
    
      case (PERIODIC_BC)
        ! Need to find which edge is closest
        ! Get position in the plane & direction
        rl = r(self % plane) - self % origin
        ul = u(self % plane)

        iHit = -1
        d = INF

        do i = 1, 6
      
          iNext = mod(i,6) + 1
    
          e = self % verts(iNext,:) - self % verts(i,:)
          b = rl - self % verts(iHit,:)

          d0 = dot_product(e, b) / dot_product(e, e)
          d0 = max(ZERO, min(ONE, d0))

          proj = self % verts(i,:) + d0 * e
          d0 = norm2(rl - proj)

          if (d0 < d) then
            d = d0
            iHit = i
          end if

        end do

        ! Find normal of the hit plane
        iNext = mod(iHit,6) + 1
        e = self % verts(iNext,:) - self % verts(iHit,:)
        n(1) = -e(2)
        n(2) = e(1)
        n = n / norm2(e)

        ! Calculate displacement and perform translation
        r(self % plane) = r(self % plane) - TWO * self % halfwidth * n

      case default
        call fatalError(Here, 'Unrecognised BC: '// numToChar(self % BC))

    end select

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
    class(hexagon), intent(in)                 :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u
    real(defReal), dimension(2)                :: rl, a1, a2, c
    integer(shortInt), dimension(2)            :: iShift
    real(defReal), dimension(2,2)              :: A, Ainv
    real(defReal)                              :: invDet
    character(100), parameter :: Here = 'transformBC (hexagon_class.f90)'

    ! Apply BC
    select case(self % BC)
      case (VACUUM_BC)
        ! Do nothing. Pass

      case (PERIODIC_BC)

        rl = r(self % plane) - self % origin

        ! Define hexagonal array vectors
        ! I THINK THESE VECTORS ARE FUCKING STUPID
        if (self % flatAxis == 1) then
          a1 = self % halfwidth * [ 3 * HALF, sqrt(3.0_defReal) * HALF]
          a2 = self % halfwidth * [ ZERO, sqrt(3.0_defReal)]
        else
          a1 = self % halfwidth * [sqrt(3.0_defReal), ZERO]
          a2 = self % halfwidth * [sqrt(3.0_defReal) * HALF, 3 * HALF]
        end if

        ! No check on determinant: magnitude of vectors can be small
        ! but it should be hard to become singular!
        invDet = ONE / (a1(1) * a2(2) - a1(2) * a2(1))

        Ainv(1,1) = a2(2) * invDet
        Ainv(1,2) = -a1(2) * invDet
        Ainv(2,1) = -a2(1) * invDet
        Ainv(2,2) = a1(1) * invDet

        ! Convert position to fractional lattice coords and round
        ! to nearest lattice shifts
        c = matmul(Ainv, rl)
        A(:,1) = a1
        A(:,2) = a2
        r(self % plane) = rl - matmul(A,nint(c))

      case default
        call fatalError(Here, 'Unrecognised BC: '// numToChar(self % BC))
    end select
    
  end subroutine transformBC

end module hexagon_class
