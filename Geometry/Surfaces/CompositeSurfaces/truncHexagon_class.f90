module truncHexagon_class

  use numPrecision
  use universalVariables
  use genericProcedures,  only : fatalError, numToChar, swap
  use dictionary_class,   only : dictionary
  use surface_inter,      only : surface, kill_super => kill

  implicit none
  private

  integer(shortInt), parameter :: pointType = 1, flatType = 2

  real(defReal), parameter :: THREE = 3.0_defReal, SQRT3 = sqrt(3.0_defReal), &
                              TWO_SQRT3 = TWO / sqrt(3.0_defReal), THREE_TWO = 1.5_defReal
  !!
  !! Finite-height hexagon aligned with x, y or z axis.
  !! The hexagon can also be oriented either to have the point or flat
  !! aligning with the first of the remaining two axis.
  !!
  !! Three different axis types are avaliable
  !!   xTruncHexagon -> parallel with X-axis
  !!   yTruncHexagon -> parallel with Y-axis
  !!   zTruncHexagon -> parallel with Z-axis
  !!
  !! One can then specify the orientation to be type 1 or 2.
  !! Type 1 has the flat surface perpendicular to the first of the two axis in the plane
  !! Type 2 has the flat surface parallel to the first of the two axis in the plane
  !!
  !! Surface Tolerance: SURF_TOL
  !!
  !! Sample Dictionary Input:
  !!   x { type xTruncHexagon; id 92; orientation 1; origin (0.0 0.0 9.0); halfwidth 7; halfheight 7;}
  !!   y { type yTruncHexagon; id 94; orientation 2; origin (0.0 0.0 9.0); halfwidth 3; halfheight 5;}
  !!
  !! Origin entry in axis along the cylinder direction is ignored.
  !! Halfwidth is half the flat-to-flat distance.
  !! Halfheight is half the top-to-bottom distance.
  !!
  !! Boundary Conditions:
  !!   Only a single BC value within the hexagon plane.
  !!   Any BC can be used along the axial direction.
  !!
  !!   For now, each radial face must have the same BC. Supports either periodic or vacuum BCs.
  !!   TODO: add reflective transform BCs.
  !!
  !! Private Members:
  !!   origin      -> position of the middle of the hexagon
  !!   halfwidth   -> Halfwidth of the hexagon from one flat to the opposite. Same in all directions
  !!   halfheight  -> Halfheight of the hexagon from the top to the bottom.
  !!   BC          -> Boundary conditions flag for all faces in the plane and on the top and bottom.
  !!   plane       -> Indices of the axis, which are the plane of the cylinder
  !!                  (e.g. for xTruncHexagon, y-z)
  !!   axis        -> Index of the axis of the plane
  !!   flatAxis    -> Index into plane, deciding which plane axis has the flat of the hexagon
  !!   pointAxis   -> Index into plane, deciding which plane axis has the point of the hexagon
  !!   verts       -> Vertices of the hexagon
  !!
  !! Interface:
  !!   surface interface
  !!
  type, public, extends(surface) :: truncHexagon
    private
    real(defReal), dimension(3)     :: origin     = ZERO
    real(defReal)                   :: halfwidth  = ZERO
    real(defReal)                   :: halfheight = ZERO
    integer(shortInt), dimension(3) :: BC = VACUUM_BC
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
    procedure :: normal
    procedure :: kill
    procedure :: setBC
    procedure :: explicitBC
    procedure :: transformBC
  end type truncHexagon


contains

  !!
  !! Return surface type name
  !!
  !! See surface_inter for more details
  !!
  pure function myType(self) result(str)
    class(truncHexagon), intent(in) :: self
    character(:), allocatable       :: str
    character(100), parameter       :: Here = 'myType (truncHexagon_class.f90)'

    select case (self % axis)
      case (X_AXIS)
        str = 'xTruncHexagon'

      case (Y_AXIS)
        str = 'yTruncHexagon'

      case (Z_AXIS)
        str = 'zTruncHexagon'

      case default
        str = 'Unknown truncated hexagon'

    end select

  end function myType

  !!
  !! Initialise truncated hexagon from a dictionary
  !!
  !! See surface_inter for more details
  !!
  !! Errors:
  !!  fatalError for -ve id or meaningful halfspace.
  !!  fatalError for unrecognised oritnetation
  !!
  subroutine init(self, dict)
    class(truncHexagon), intent(inout)       :: self
    class(dictionary), intent(in)            :: dict
    integer(shortInt)                        :: id, N, i, orient
    real(defReal), dimension(:), allocatable :: temp
    character(nameLen)                       :: type
    real(defReal), dimension(6)              :: angles
    character(100), parameter :: Here = 'init (truncHexagon_class.f90)'

    ! Load id
    call dict % get(id,'id')
    if (id <= 0) call fatalError(Here, 'ID must be <= 0. Is: '//numToChar(id))
    call self % setID(id)

    ! Select type
    call dict % get(type,'type')
    select case(type)
      case('xTruncHexagon')
        self % axis = X_AXIS
        self % plane = [Y_AXIS, Z_AXIS]

      case('yTruncHexagon')
        self % axis = Y_AXIS
        self % plane = [X_AXIS, Z_AXIS]

      case('zTruncHexagon')
        self % axis = Z_AXIS
        self % plane = [X_AXIS, Y_AXIS]

      case default
        call fatalError(Here, 'Unknown type of truncated hexagon: '//type)

    end select
    
    ! Load origin
    call dict % get(temp,'origin')
    N = size(temp)
    if (N /= 3) call fatalError(Here,'origin must have size 3. Has: '//numToChar(N))
    self % origin = temp

    ! Load halfwidth
    call dict % get(self % halfwidth,'halfwidth')
    if (self % halfwidth <= ZERO) call fatalError(Here, 'halfwidth cannot have a zero/-ve value.')

    ! Load halfheight
    call dict % get(self % halfheight,'halfheight')
    if (self % halfheight <= ZERO) call fatalError(Here, 'halfheight cannot have a zero/-ve value.')

    ! Get orientation
    call dict % get(orient, 'orientation')
    if (all(orient /= [pointType, flatType])) then
      call fatalError(Here, 'Unrecognised hexagon orientation: '//numToChar(orient))
    end if

    if (orient == pointType) then
      self % flatAxis = 1
      self % pointAxis = 2
      angles = [(i, i = 0,5)] * PI/3 + PI/6
    else
      self % flatAxis = 2
      self % pointAxis = 1
      angles = [(i, i = 0,5)] * PI/3
    end if
    self % verts(:,1) = TWO_SQRT3 * self % halfwidth * cos(angles)
    self % verts(:,2) = TWO_SQRT3 * self % halfwidth * sin(angles)

  end subroutine init

  !!
  !! Return axis-aligned bounding box for the surface
  !!
  !! See surface_inter for details
  !!
  pure function boundingBox(self) result(aabb)
    class(truncHexagon), intent(in)  :: self
    real(defReal), dimension(6)      :: aabb
    real(defReal)                    :: originF, originP, originA
    integer(shortInt)                :: f, p

    ! Two of the planes will be parallel to either the first or second axis in the plane
    f = self % plane(self % flatAxis)
    p = self % plane(self % pointAxis)
    originF = self % origin(f)
    originP = self % origin(p)
    originA = self % origin(self % axis)

    aabb(f)     = originF - self % halfwidth
    aabb(f + 3) = originF + self % halfwidth
      
    ! The other two are shifted further by 2/sqrt(3)
    aabb(p)     = originP - self % halfwidth * TWO_SQRT3
    aabb(p + 3) = originP + self % halfwidth * TWO_SQRT3

    aabb(self % axis)      = originA - self % halfheight
    aabb(self % axis + 3)  = originA + self % halfheight

  end function boundingBox

  !!
  !! Evaluate surface expression c = F(r)
  !!
  !! See surface_inter for details
  !!
  pure function evaluate(self, r) result(c)
    class(truncHexagon), intent(in)         :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: c
    real(defReal), dimension(3)             :: rb
    real(defReal)                           :: p, q, z

    ! Move to origin-frame and evaluate
    rb = abs(r - self % origin)
    
    p = rb(self % plane(self % flatAxis))
    q = rb(self % plane(self % pointAxis))
    z = rb(self % axis)

    c = max(p - self % halfwidth, &
            (p + SQRT3 * q) * HALF - self % halfwidth, &
            z - self % halfheight)

  end function evaluate

  !!
  !! Return distance to the surface
  !!
  !! See surface_inter for details
  !!
  pure function distance(self, r, u) result(d)
    class(truncHexagon), intent(in)         :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    real(defReal)                           :: d
    real(defReal), dimension(2)             :: rl, ul, a1, a2, e, n
    real(defReal)                           :: far, near, uProj, rProj, rz, &
                                               a_near, a_far, test_near, test_far
    integer(shortInt)                       :: i
    real(defReal), parameter :: FP_MISS_TOL = ONE + 10.0_defReal * epsilon(ONE)

    ! Get position and direction in the plane and centre
    rl = r(self % plane) - self % origin(self % plane)
    ul = u(self % plane)

    ! Initialise intersection set
    far  = huge(far)
    near = -huge(near)

    ! Loop over the 3 sets of parallel planes
    do i = 1, 3
      ! Get the normal for this pair of planes
      a1 = self % verts(i, :)
      a2 = self % verts(i + 1, :)
      e = a2 - a1
      n = [e(2), -e(1)]
      n = n / norm2(n)

      ! Project onto the normal
      uProj = dot_product(ul,n)
      rProj = dot_product(rl,n)

      a_far = sign(self % halfwidth, uProj)
      a_near = -a_far
    
      if (uProj /= ZERO ) then
        test_near = (a_near - rProj) / uProj
        test_far = (a_far - rProj) / uProj

      else ! Line is either fully inside or fully outside
        test_near = sign(INF, (a_near - rProj))
        test_far  = sign(INF, (a_far - rProj))

        ! If direction is -0.0 (intead of just 0.0) Order may be
        ! wrong and it is necesary to swap
        if (test_near > test_far) call swap(test_near, test_far)

      end if
      far = min(far, test_far)
      near = max(near, test_near)

    end do

    ! Also check axial distances
    rz = r(self % axis) - self % origin(self % axis)
    if (u(self % axis) /= ZERO) then
      test_near = (-self % halfheight - rz) / u(self % axis)
      test_far = (self % halfheight - rz) / u(self % axis)

    else
      test_near = sign(INF, (-self % halfheight - rz))
      test_far = sign(INF, (self % halfheight - rz))

    end if
    
    ! Ensure correct order for any orientation
    if (test_far < test_near) call swap(test_far, test_near)

    ! Get intersection
    far = min(far, test_far)
    near = max(near, test_near)

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
  !!   For parallel direction halfspace is assigned by the sign of `evaluate` result.
  !!
  pure function going(self, r, u) result(halfspace)
    class(truncHexagon), intent(in)              :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    logical(defBool)                        :: halfspace
    real(defReal), dimension(2)             :: rl, ul, e, n, nTrial
    real(defReal)                           :: d, dist, proj, rz
    integer(shortInt)                       :: i

    ! Get position in the plane & direction
    rl = r(self % plane) - self % origin(self % plane)
    ul = u(self % plane)

    ! Which set of planes is being crossed?
    d = -ONE
    do i = 1,3
    
      ! Get plane normal
      e = self % verts(i+1, :) - self % verts(i, :)
      nTrial = [e(2), -e(1)]
      nTrial = nTrial /norm2(nTrial)
      dist = abs(dot_product(rl,nTrial)) - self % halfwidth

      if (dist > d) then
        d = dist
        n = nTrial
      end if

    end do

    n = n * sign(ONE, rl*n)
    proj = dot_product(ul, n)

    ! Also check axial direction
    rz = abs(r(self % axis) - self % origin(self % axis))
    dist = rz - self % halfheight
    if (dist > d) then
      proj = rz * u(self % axis)
    end if

    halfspace = proj > ZERO

    ! Parallel direction
    ! Need to use position to determine halfspace
    if (proj == ZERO) then
      halfspace = self % evaluate(r) >= ZERO
    end if

  end function going
  
  !!
  !! Return normal of the surface closest to the particle
  !!
  !! See surface_inter for details
  !!
  pure function normal(self, r, u) result(n)
    class(truncHexagon), intent(in)         :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    real(defReal), dimension(3)             :: n
    real(defReal), dimension(2)             :: rl, e, nTrial, nBest
    real(defReal)                           :: d, dist, rz
    integer(shortInt)                       :: i
    real(defReal), dimension(3)             :: n0
    
    ! Get position in the plane & direction
    rl = r(self % plane) - self % origin(self % plane)

    ! Which set of planes is being crossed?
    d = -INF
    do i = 1,3
    
      ! Get plane normal
      e = self % verts(i+1, :) - self % verts(i, :)
      nTrial = [e(2), -e(1)]
      nTrial = nTrial /norm2(nTrial)
      dist = abs(dot_product(rl,nTrial)) - self % halfwidth

      if (dist > d) then
        d = dist
        nBest = nTrial

      ! Catch corners
      elseif (dist == d) then
        nBest = nBest + nTrial
        nBest = nBest / norm2(nBest)
      end if

    end do

    n(self % plane) = nBest * sign(ONE, rl * nBest)
    n(self % axis) = ZERO

    ! Also check axial planes
    rz = r(self % axis) - self % origin(self % axis)
    if ((abs(rz) - self % halfheight) > d) then
      n = ZERO
      n(self % axis) = sign(ONE, rz)
    elseif (abs(rz) - self % halfheight == d) then
      ! Catch corners
      n0 = ZERO
      n0(self % axis) = sign(ONE, rz)
      n = n + n0
      n = n /norm2(n)
    end if

  end function normal

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(truncHexagon), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Local
    self % origin      = ZERO
    self % halfwidth   = ZERO
    self % halfheight  = ZERO
    self % BC          = VACUUM_BC
    self % plane       = 0
    self % flatAxis    = 0
    self % pointAxis   = 0
    self % axis        = 0
    self % verts       = ZERO

  end subroutine kill

  !!
  !! Set boundary conditions
  !!
  !! See surface_inter for details
  !!
  subroutine setBC(self, BC)
    class(truncHexagon), intent(inout)          :: self
    integer(shortInt), dimension(:), intent(in) :: BC
    integer(shortInt)                           :: i
    character(100),parameter                    :: Here = 'setBC (hexagon_class.f90)'

    if(size(BC) /= 6) call fatalError(Here, 'Wrong size of BC string. Must be 6.')
    if(any(BC(self % plane) /= BC(self % plane + 3))) call fatalError(Here, 'BCs in the plane must be the same.')
    if((BC(self % axis) == PERIODIC_BC .and. BC(self % axis + 3) /= PERIODIC_BC) .or. &
       (BC(self % axis) /= PERIODIC_BC .and. BC(self % axis + 3) == PERIODIC_BC)) then
      call fatalError(Here,'Must have matching periodic boundary conditions in the axial direction')
    end if

    ! Load BC codes
    self % BC(1) = BC(self % plane(1))
    self % BC(2) = BC(self % axis)
    self % BC(3) = BC(self % axis+3)

    ! Verify that all BC flags make sense
    select case(self % BC(1))
      case (VACUUM_BC, PERIODIC_BC)
          ! Do nothing, pass

      ! Reflective not fully supported
      ! TODO: Implement transform reflective BCs
      case (REFLECTIVE_BC)
        call fatalError(Here,'Reflective boundaries not supported in the plane.')
      case default
        call fatalError(Here,'Unrecognised BC: '//numToChar(BC))

    end select
    
    do i = 2,3
      select case(self % BC(i))
        case (VACUUM_BC, PERIODIC_BC, REFLECTIVE_BC)
            ! Do nothing, pass

        case default
          call fatalError(Here,'Unrecognised BC: '//numToChar(BC))

      end select
    end do

  end subroutine setBC

  !!
  !! Apply explicit BCs
  !!
  !! See surface_inter for details
  !!
  !! Note:
  !!   - Go through all directions in order to account for corners
  !!   - The periodic shift in a corner is decided by the particle direction
  !!
  subroutine explicitBC(self, r, u)
    class(truncHexagon), intent(in)            :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u
    integer(shortInt)                          :: i, iPlane, bc
    real(defReal)                              :: dist, proj, maxProj, rz
    real(defReal), dimension(2)                :: rl, e, v1, v2, n, nBest
    character(100), parameter :: Here = 'explicitBC (truncHexagon_class.f90)'

    ! Get position in the plane & direction
    rl = r(self % plane) - self % origin(self % plane)

    maxProj = -INF
    iPlane = 0

    ! Identify which face is being crossed by the particle, including direction
    ! This tie-break is necessary when crossing corners.
    ! Otherwise, shifting would occur based on loop-ordering.
    ! For a particle pointing in between the two normals, the shift
    ! will indeed be determined by loop-ordering.
    do i = 1, 3
      v1 = self % verts(i, :)
      v2 = self % verts(i + 1, :)
      e = v2 - v1
      n = [e(2), -e(1)]
      n = n / norm2(n)

      dist = dot_product(rl, n)

      ! Are we on the plane?
      if (abs(dist) >= self % halfwidth - self % surfTol()) then
        
        n = n * sign(ONE, dist)
        proj = dot_product(u(self % plane), n)
        
        ! The face with the largest positive projection is the true exit face
        if (proj > maxProj) then
          maxProj = proj
          iPlane = i
          nBest = n
        end if
      end if
    end do

    ! Check the other axis
    rz = r(self % axis) - self % origin(self % axis)
    if (abs(rz) >= self % halfheight - self % surfTol()) then

      proj = u(self % axis) * sign(ONE, rz)

      if (proj > maxProj) then
        maxProj = proj
        iPlane = -1
      end if
    end if

    ! Apply the BC only once, using the physically correct face
    ! Hits the hexagon planes
    if (iPlane > 0) then
      select case(self % BC(1))
        case (PERIODIC_BC)
          r(self % plane) = r(self % plane) - TWO * self % halfwidth * nBest

        ! This should never be called at present!
        ! Need to implement transform BCs!
        case (REFLECTIVE_BC)
          u(self % plane) = u(self % plane) - TWO * dot_product(nBest, u(self % plane)) * nBest

        case (VACUUM_BC)
          ! Do nothing

        case default
          call fatalError(Here, 'Unrecognised BC: '// numToChar(self % BC(1)))
      end select

    ! Hits the top/bottom planes
    elseif (iPlane == -1) then
      ! Top or bottom?
      if (rz > 0) then
        bc = self % BC(3)
      else
        bc = self % BC(2)
      end if

      select case(bc)
        case (PERIODIC_BC)
          r(self % axis) = r(self % axis) - TWO * self % halfheight * sign(ONE, rz)

        case (REFLECTIVE_BC)
          u(self % axis) = u(self % axis) - TWO * u(self % axis)

        case (VACUUM_BC)
          ! Do nothing

        case default
          call fatalError(Here, 'Unrecognised BC: '// numToChar(bc))
      end select

    end if

    ! Otherwise, no surface was hit and BCs are not applied

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
    class(truncHexagon), intent(in)            :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u
    real(defReal), dimension(2)                :: rl
    integer(shortInt), dimension(2)            :: latShift
    real(defReal), dimension(2,2)              :: Ainv
    real(defReal)                              :: invDet, a_bar, a0, d, r0
    real(defReal), dimension(2)                :: latCoord, a1, a2
    integer(shortInt)                          :: Ri, t, bc, ax
    character(100), parameter :: Here = 'transformBC (truncHexagon_class.f90)'

    ! Apply BC in the hexagon plane
    select case(self % BC(1))
      case (VACUUM_BC)
        ! Do nothing. Pass

      case (PERIODIC_BC)

        rl = r(self % plane) - self % origin(self % plane)
        ! Define the lattice basis vectors 
        ! For a flat-to-flat distance of 2H, the distance between centers 
        ! of adjacent hexagons is SQRT(3) * H along the flat-axis.
        if (self % flatAxis == 1) then
          a1 = [ TWO * self % halfwidth, ZERO ]
          a2 = [ self % halfwidth, SQRT3 * self % halfwidth ]
        else
          a1 = [ SQRT3 * self % halfwidth, self % halfwidth ]
          a2 = [ ZERO, TWO * self % halfwidth ]
        end if
        
        ! No check on determinant: magnitude of vectors can be small
        ! but it should be hard to become singular provided halfwidth is positive!
        invDet = ONE / (a1(1) * a2(2) - a1(2) * a2(1))

        Ainv(1,1) = a2(2) * invDet
        Ainv(1,2) = -a2(1) * invDet
        Ainv(2,1) = -a1(2) * invDet
        Ainv(2,2) = a1(1) * invDet

        ! Convert to lattice coordinates
        latCoord = matmul(Ainv, rl)
        latShift = nint(latCoord)

        rl = rl - (real(latShift(1), defReal) * a1 + real(latShift(2), defReal) * a2)

        r(self % plane) = rl + self % origin(self % plane)
      
      ! This should never be called at present!
      case (REFLECTIVE_BC)
        call fatalError(Here, 'Reflective boundary conditions not supported.')
      case default
        call fatalError(Here, 'Unrecognised BC: '// numToChar(self % BC))
    end select

    ! Apply BC in the axial direction
    ! Calculate halfwidth reduced by the surface_tolerance
    ! Necessary to capture particles at the boundary
    a_bar = self % halfheight - self % surfTol()

    ! Calculate distance (in # of transformations) in axial direction
    ax = self % axis
    Ri = ceiling(abs(r(self % axis) - self % origin(ax)) / a_bar) / 2

    trans : do t = 1, Ri
      ! Find position wrt origin
      r0 = r(ax) - self % origin(ax)

      ! Choose correct BC
      if (r0 < ZERO) then
        bc = self % BC(2)
      else
        bc = self % BC(3)
      end if

      ! Apply BC
      select case(bc)
        case (VACUUM_BC)
          ! Do nothing. Pass

        case (REFLECTIVE_BC)
          ! Find position of the plane
          a0 = sign(self % halfheight, r0) + self % origin(ax)

          ! Calculate displacment and perform reflection
          d = r(ax) - a0
          r(ax) = r(ax) - TWO * d
          u(ax) = -u(ax)

        case (PERIODIC_BC)
          ! Calculate displacement and perform translation
          d = sign(self % halfheight, r0)
          r(ax) = r(ax) - TWO * d

        case default
          call fatalError(Here, 'Unrecognised BC: '// numToChar(bc))
      end select
    end do trans
    
  end subroutine transformBC

end module truncHexagon_class
