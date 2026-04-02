module hexagon_class

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
  !! Hexagon aligned with x, y or z axis.
  !! The hexagon can also be oriented either to have the point or flat
  !! aligning with the first of the remaining two axis.
  !!
  !! Three different axis types are avaliable
  !!   xHexagon -> parallel with X-axis
  !!   yHexagon -> parallel with Y-axis
  !!   zHexagon -> parallel with Z-axis
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
  !! Origin entry in axis along the cylinder direction is ignored.
  !! Halfwidth is the flat-to-flat distance.
  !!
  !! Boundary Conditions:
  !!   Only a single BC value.
  !!
  !!   For now, each face must have the same BC. Supports either periodic or vacuum BCs.
  !!   TODO: add reflective transform BCs.
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
  !!   verts       -> Vertices of the hexagon
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
    procedure :: normal
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
        call fatalError(Here, 'Unknown type of hexagon: '//type)

    end select
    
    ! Load origin
    call dict % get(temp,'origin')
    N = size(temp)
    if (N /= 3) call fatalError(Here,'origin must have size 3. Has: '//numToChar(N))
    self % origin = temp(self % plane)

    ! Load halfwidth
    call dict % get(self % halfwidth,'halfwidth')
    if (self % halfwidth <= ZERO) call fatalError(Here, 'halfwidth cannot have a zero/-ve value.')

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
    real(defReal)                           :: p, q

    ! Move to origin-frame and evaluate
    rb = abs(r(self % plane) - self % origin)
    
    p = rb(self % flatAxis)
    q = rb(self % pointAxis)

    c = max(p - self % halfwidth, &
            (p + SQRT3 * q) * HALF - self % halfwidth)

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
    real(defReal), dimension(2)             :: rl, ul, a1, a2, e, n
    real(defReal)                           :: far, near, uProj, rProj, &
                                               a_near, a_far, test_near, test_far
    integer(shortInt)                       :: i
    real(defReal), parameter :: FP_MISS_TOL = ONE + 10.0_defReal * epsilon(ONE)

    ! Get position and direction in the plane and centre
    rl = r(self % plane) - self % origin
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
    class(hexagon), intent(in)              :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    logical(defBool)                        :: halfspace
    real(defReal), dimension(2)             :: rl, ul, e, n, nTrial
    real(defReal)                           :: d, dist, proj
    integer(shortInt)                       :: i

    ! Get position in the plane & direction
    rl = r(self % plane) - self % origin
    ul = u(self % plane)

    ! Which set of planes is being crossed?
    d = -ONE
    do i = 1,3
    
      ! Get plane normal
      e = self % verts(i+1, :) - self % verts(i, :)
      nTrial = [e(2), -e(1)]
      nTrial = nTrial /norm2(nTrial)
      dist = abs(dot_product(rl,nTrial))

      if (dist > d) then
        d = dist
        n = nTrial
      end if

    end do

    n = n * sign(ONE, rl*n)
    proj = dot_product(ul, n)

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
    class(hexagon), intent(in)              :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    real(defReal), dimension(3)             :: n
    real(defReal), dimension(2)             :: rl, e, nTrial, nBest
    real(defReal)                           :: d, dist
    integer(shortInt)                       :: i
    
    ! Get position in the plane & direction
    rl = r(self % plane) - self % origin

    ! Which set of planes is being crossed?
    d = -ONE
    do i = 1,3
    
      ! Get plane normal
      e = self % verts(i+1, :) - self % verts(i, :)
      nTrial = [e(2), -e(1)]
      nTrial = nTrial /norm2(nTrial)
      dist = abs(dot_product(rl,nTrial))

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

  end function normal

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
    character(100),parameter                    :: Here = 'setBC (hexagon_class.f90)'

    if(size(BC) /= 6) call fatalError(Here, 'Wrong size of BC string. Must be 6')
    if(all(BC /= BC(1))) call fatalError(Here, 'All BCs must be identical')

    ! Load BC codes
    self % BC = BC(1)

    ! Verify that all BC flags make sense
    select case(self % BC)
      case (VACUUM_BC, PERIODIC_BC)
          ! Do nothing, pass

      ! Reflective not fully supported
      ! TODO: Implement transform reflective BCs
      case (REFLECTIVE_BC)
        call fatalError(Here,'Reflective boundaries not supported.')
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
  !!   - The periodic shift in a corner is decided by the particle direction
  !!
  subroutine explicitBC(self, r, u)
    class(hexagon), intent(in)                 :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u
    integer(shortInt)                          :: i, iPlane
    real(defReal)                              :: dist, proj, maxProj
    real(defReal), dimension(2)                :: rl, e, v1, v2, n, nBest
    character(100), parameter :: Here = 'explicitBC (hexagon_class.f90)'

    ! Get position in the plane & direction
    rl = r(self % plane) - self % origin

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

    ! Apply the BC only once, using the physically correct face
    if (iPlane > 0) then
      select case(self % BC)
        case (PERIODIC_BC)
          r(self % plane) = r(self % plane) - TWO * self % halfwidth * nBest

        ! This should never be called at present!
        ! Need to implement transform BCs!
        case (REFLECTIVE_BC)
          u(self % plane) = u(self % plane) - TWO * dot_product(nBest, u(self % plane)) * nBest

        case (VACUUM_BC)
          ! Do nothing

        case default
          call fatalError(Here, 'Unrecognised BC: '// numToChar(self % BC))
      end select
    end if

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
    real(defReal), dimension(2)                :: rl
    integer(shortInt), dimension(2)            :: baseShift, latShift
    real(defReal), dimension(2,2)              :: Ainv
    real(defReal)                              :: invDet, dist2, minDist2
    real(defReal), dimension(2)                :: latCoord, a1, a2, centre
    integer(shortInt)                          :: i, j
    character(100), parameter :: Here = 'transformBC (hexagon_class.f90)'

    ! Apply BC
    select case(self % BC)
      case (VACUUM_BC)
        ! Do nothing. Pass

      case (PERIODIC_BC)

        rl = r(self % plane) - self % origin
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
        baseShift = nint(latCoord)

        ! We check a 3x3 local grid of lattice points around our guess
        ! to see which center is strictly closest in Cartesian space.
        minDist2 = INF
        latShift = baseShift

        do i = baseShift(1) - 1, baseShift(1) + 1
          do j = baseShift(2) - 1, baseShift(2) + 1

            ! Calculate Cartesian position of this neighboring lattice center
            centre = i * a1 + j * a2

            ! Squared distance from particle to this lattice center
            dist2 = sum((rl - centre)**2)

            if (dist2 < minDist2) then
              minDist2 = dist2
              latShift = [i, j]
            end if

          end do
        end do

        rl = rl - latShift(1) * a1 - latShift(2) * a2

        r(self % plane) = rl + self % origin
      
      ! This should never be called at present!
      case (REFLECTIVE_BC)
        call fatalError(Here, 'Reflective boundary conditions not supported.')
      case default
        call fatalError(Here, 'Unrecognised BC: '// numToChar(self % BC))
    end select
    
  end subroutine transformBC

end module hexagon_class
