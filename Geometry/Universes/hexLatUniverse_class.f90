module hexLatUniverse_class

  use numPrecision
  use universalVariables, only : INF, SURF_TOL, X_AXIS, Y_AXIS, Z_AXIS
  use genericProcedures,  only : fatalError, numToChar, swap
  use dictionary_class,   only : dictionary
  use coord_class,        only : coord
  use charMap_class,      only : charMap
  use surfaceShelf_class, only : surfaceShelf
  use box_class,          only : box
  use cell_inter,         only : cell
  use cellShelf_class,    only : cellShelf
  use universe_inter,     only : universe, kill_super => kill, charToFill

  implicit none
  private

  ! Orientation flags (same convention as hexagon_class/truncHexagon_class)
  integer(shortInt), parameter :: pointType = 1, flatType = 2

  ! Options for the offset map (same convention as latUniverse_class)
  integer(shortInt), parameter :: local = 1, noOffset = 0

  real(defReal), parameter :: SQRT3 = sqrt(3.0_defReal), TWO_SQRT3 = TWO / sqrt(3.0_defReal)

  ! Local surface indexes
  ! Note: FACE<n>_POS/NEG refer to the 3 pairs of parallel edges of a hexagonal cell.
  ! Public so they can be accessed by unit tests
  integer(shortInt), public, parameter :: FACE1_POS = -1, FACE1_NEG = -2, &
                                          FACE2_POS = -3, FACE2_NEG = -4, &
                                          FACE3_POS = -5, FACE3_NEG = -6, &
                                          AX_MIN = -7, AX_MAX = -8, OUTLINE_SURF = -9

  !!
  !! Hexagonal lattice universe, stacked in layers along a Cartesian axis
  !!
  !! Universe consists of a hexagonal arrangement of cells (e.g. 5 x 5), which may be
  !! further stacked into layers along the axis perpendicular to the hexagonal plane
  !! (e.g. giving a 5 x 5 x 2 lattice). Centre of the lattice is placed at the origin.
  !! An additional cell is placed beyond the lattice called background (or out) cell.
  !!
  !! The lattice is defined using oblique (axial) co-ordinates i,j in the hexagonal plane.
  !! Cell centres are placed at: r_ij = (i - (N1+1)/2) * a1 + (j - (N2+1)/2) * a2
  !! where a1 & a2 are the two lattice basis vectors (60 deg apart, magnitude equal to pitch),
  !! which depend on the chosen orientation of the individual hexagonal cells. This gives a
  !! parallelogram-shaped map of hexagons, in the style of a Serpent x-type/y-type hexagonal
  !! lattice.
  !!
  !! Two orientations are available (matching hexagon_class):
  !!   1 -> point,  the hexagon has a vertex aligned with the 1st in-plane axis
  !!   2 -> flat,   the hexagon has a flat face aligned with the 1st in-plane axis
  !!
  !! Local ID is 1 at the (i,j,k) = (1,1,1) cell. It increases first with i, then j and
  !! lastly k (the axial layer). Cells inside the lattice can only be filled with a universe
  !! (given as integer ID). Background cell can have any filling given by keyword (material or
  !! universe).
  !!
  !! Every lattice cell has an offset to its centre (so the centre of the nested universe is in
  !! the centre of the lattice cell). Optionally an offset map can be provided, analogous to the
  !! one available in latUniverse, determining whether to apply an offset in a given cell
  !! position.
  !!
  !! Minimum pitch is set to 10 * SURF_TOL
  !!
  !! Sample Input Dictionary (3D):
  !!   hlatt { id 7;
  !!           type zHexLattice;
  !!           #origin (0.0 0.0 0.0); #
  !!           #rotation (0.0 0.0 0.0); #
  !!           orientation 1;
  !!           pitch (1.5 3.0);       // (radial pitch, axial pitch)
  !!           shape (3 3 2);         // (N1, N2, N3)
  !!           padMat void;
  !!           map ( 1 2 3            // Top layer, back row
  !!                 4 5 6            // Top layer, middle row
  !!                 7 8 9            // Top layer, front row
  !!
  !!                10 11 12
  !!                13 14 15
  !!                16 17 18  );
  !!   }
  !!
  !! Sample Input Dictionary (2D):
  !!   hlatt2D { id 8;
  !!             type zHexLattice;
  !!             orientation 2;
  !!             shape (2 2 0);       // 0 indicates infinite extent along the axis
  !!             pitch (1.0 0.0);     // Axial pitch is ignored, use 0.0 for clarity
  !!             padMat void;
  !!             map ( 1 2
  !!                   2 1);
  !!   }
  !!
  !! NOTE: Input in MAP for a single layer is WYSIWYG, as in latUniverse: the first row of the
  !!   map is the row with the highest j (2nd in-plane axial co-ordinate).
  !!
  !! Private Members:
  !!   axis       -> Index of the axis perpendicular to the hexagonal plane
  !!   plane      -> Indices of the 2 axes forming the hexagonal plane
  !!   flatAxis   -> Index into plane, deciding which plane axis has the flat of the hexagon
  !!   pointAxis  -> Index into plane, deciding which plane axis has the point of the hexagon
  !!   halfwidth  -> Halfwidth (flat-to-flat) of an individual hexagonal cell
  !!   axialPitch -> Pitch of the layers along the axis
  !!   sizeN      -> Number of cells: (N1, N2, N3) with N1 & N2 in-plane and N3 axial
  !!   a1, a2     -> Lattice basis vectors of the hexagonal tiling (in-plane, 2D)
  !!   Ainv       -> Inverse of the matrix [a1 a2], used to convert Cartesian to oblique co-ords
  !!   faceNormal -> Unit normals of the 3 pairs of parallel edges of a hexagonal cell
  !!   faceOffset -> For each of the 3 face pairs, the (di,dj) neighbour shift on the +normal side
  !!   axialCorner-> Position of the minimum edge of the lowest axial layer
  !!   outline    -> Box type surface, bounding (from outside) the valid extent of the lattice
  !!   outLocalID -> LocalID of the background cell
  !!   offset     -> Flag to disable all offsets
  !!   offsetMap  -> Map determining which cells have a lattice offset or not
  !!
  !! Interface:
  !!   universe interface
  !!
  type, public, extends(universe) :: hexLatUniverse
    private
    integer(shortInt)                 :: axis      = 0
    integer(shortInt), dimension(2)   :: plane     = 0
    integer(shortInt)                 :: flatAxis  = 0
    integer(shortInt)                 :: pointAxis = 0
    real(defReal)                     :: halfwidth  = ZERO
    real(defReal)                     :: axialPitch = ZERO
    integer(shortInt), dimension(3)   :: sizeN      = 0
    real(defReal), dimension(2)       :: a1 = ZERO, a2 = ZERO
    real(defReal), dimension(2,2)     :: Ainv       = ZERO
    real(defReal), dimension(3,2)     :: faceNormal = ZERO
    integer(shortInt), dimension(3,2) :: faceOffset = 0
    real(defReal)                     :: axialCorner = ZERO
    type(box)                         :: outline
    integer(shortInt)                 :: outLocalID = 0
    logical(defBool)                  :: offset = .true.
    integer(shortInt), dimension(:), allocatable :: offsetMap
  contains
    ! Superclass procedures
    procedure :: init
    procedure :: kill
    procedure :: findCell
    procedure :: distance
    procedure :: cross
    procedure :: cellOffset
    procedure :: getNormal

    ! Private helpers
    procedure, private :: cellCentre
    procedure, private :: nearestIJ
  end type hexLatUniverse

contains

  !!
  !! Initialise Universe
  !!
  !! See universe_inter for details.
  !!
  !! Errors:
  !!   fatalError if input is invalid.
  !!
  subroutine init(self, fill, dict, cells, surfs, mats)
    class(hexLatUniverse), intent(inout)                      :: self
    integer(shortInt), dimension(:), allocatable, intent(out) :: fill
    class(dictionary), intent(in)                             :: dict
    type(cellShelf), intent(inout)                            :: cells
    type(surfaceShelf), intent(inout)                         :: surfs
    type(charMap), intent(in)                                 :: mats
    real(defReal), dimension(:), allocatable       :: temp
    integer(shortInt), dimension(:), allocatable   :: tempI
    integer(shortInt)                              :: N, i, j, outFill, val, orient, f, d, bestIdx
    type(dictionary)                                :: tempDict
    integer(shortInt), dimension(:,:), allocatable :: tempMap
    character(nameLen)                             :: name, type
    real(defReal), dimension(6)                    :: angles
    real(defReal), dimension(6,2)                  :: verts
    real(defReal), dimension(2)                    :: edge, edgeN
    real(defReal), dimension(3,2)                  :: dirs
    integer(shortInt), dimension(3,2)              :: offs
    real(defReal)                                  :: invDet, dotv, bestDot, hw
    real(defReal), dimension(2)                    :: c1, c2, c3, c4
    real(defReal), dimension(3)                    :: halfwidth3
    character(100), parameter :: Here = 'init (hexLatUniverse_class.f90)'

    ! Setup the base class
    ! With: id, origin, rotation...
    call self % setupBase(dict)

    ! Perform offsets on every cell?
    call dict % getOrDefault(self % offset, 'offset', .true.)

    ! Select type -> axis & plane
    call dict % get(type, 'type')
    select case(type)
      case('xHexLattice')
        self % axis  = X_AXIS
        self % plane = [Y_AXIS, Z_AXIS]

      case('yHexLattice')
        self % axis  = Y_AXIS
        self % plane = [X_AXIS, Z_AXIS]

      case('zHexLattice')
        self % axis  = Z_AXIS
        self % plane = [X_AXIS, Y_AXIS]

      case default
        call fatalError(Here, 'Unknown type of hexLatUniverse: '//type)

    end select

    ! Load pitch: (radial flat-to-flat pitch, axial pitch)
    call dict % get(temp, 'pitch')
    N = size(temp)
    if (N /= 2) call fatalError(Here, 'Pitch must have size 2 (radial, axial). Has: '//numToChar(N))
    if (temp(1) <= ZERO) call fatalError(Here, 'Radial pitch cannot have a zero/-ve value.')
    self % halfwidth  = temp(1) * HALF
    self % axialPitch = temp(2)

    ! Load Shape: (N1, N2, N3)
    call dict % get(tempI, 'shape')
    N = size(tempI)
    if (N /= 3) call fatalError(Here, 'Shape must have size 3. Has: '//numToChar(N))
    if (any(tempI(1:2) < 1)) call fatalError(Here, 'In-plane shape entries must be >= 1.')
    if (tempI(3) < 0) call fatalError(Here, 'Axial shape entry cannot be -ve.')
    self % sizeN = tempI

    ! Detect reduced (infinite) axial dimension
    if (self % sizeN(3) == 0) then
      self % sizeN(3)   = 1
      self % axialPitch = TWO * INF
    end if

    ! Check for invalid pitch
    if (self % halfwidth * TWO < 10 * SURF_TOL) then
      call fatalError(Here, 'Radial pitch must be larger than: '//numToChar(10 * SURF_TOL))
    else if (self % axialPitch < 10 * SURF_TOL) then
      call fatalError(Here, 'Axial pitch must be larger than: '//numToChar(10 * SURF_TOL))
    end if

    ! Get orientation
    call dict % get(orient, 'orientation')
    if (all(orient /= [pointType, flatType])) then
      call fatalError(Here, 'Unrecognised hexagon orientation: '//numToChar(orient))
    end if

    hw = self % halfwidth
    if (orient == pointType) then
      self % flatAxis  = 1
      self % pointAxis = 2
      angles = [(i, i = 0,5)] * PI/3 + PI/6
      self % a1 = [TWO * hw, ZERO]
      self % a2 = [hw, SQRT3 * hw]
    else
      self % flatAxis  = 2
      self % pointAxis = 1
      angles = [(i, i = 0,5)] * PI/3
      self % a1 = [SQRT3 * hw, hw]
      self % a2 = [ZERO, TWO * hw]
    end if
    verts(:,1) = TWO_SQRT3 * hw * cos(angles)
    verts(:,2) = TWO_SQRT3 * hw * sin(angles)

    ! Inverse of the lattice basis matrix [a1 a2]
    invDet = ONE / (self % a1(1) * self % a2(2) - self % a1(2) * self % a2(1))
    self % Ainv(1,1) =  self % a2(2) * invDet
    self % Ainv(1,2) = -self % a2(1) * invDet
    self % Ainv(2,1) = -self % a1(2) * invDet
    self % Ainv(2,2) =  self % a1(1) * invDet

    ! Determine, for each of the 3 pairs of parallel hexagon edges, the unit normal and the
    ! (i,j) neighbour shift associated with crossing that edge on its +ve normal side.
    dirs(1,:) = self % a1 / norm2(self % a1)
    dirs(2,:) = self % a2 / norm2(self % a2)
    dirs(3,:) = (self % a1 - self % a2) / norm2(self % a1 - self % a2)
    offs(1,:) = [1, 0]
    offs(2,:) = [0, 1]
    offs(3,:) = [1, -1]

    do f = 1, 3
      edge  = verts(f+1,:) - verts(f,:)
      edgeN = [edge(2), -edge(1)]
      edgeN = edgeN / norm2(edgeN)

      bestDot = ZERO
      bestIdx = 1
      do d = 1, 3
        dotv = dot_product(edgeN, dirs(d,:))
        if (abs(dotv) > abs(bestDot)) then
          bestDot = dotv
          bestIdx = d
        end if
      end do

      self % faceNormal(f,:) = edgeN
      if (bestDot > ZERO) then
        self % faceOffset(f,:) = offs(bestIdx,:)
      else
        self % faceOffset(f,:) = -offs(bestIdx,:)
      end if
    end do

    ! Axial centring offset (position of the minimum edge of the lattice)
    self % axialCorner = -(self % sizeN(3) * HALF * self % axialPitch)

    ! Calculate local ID of the background
    self % outLocalID = product(self % sizeN) + 1

    ! Build a safe (superset) bounding-box outline for the background region.
    ! The 4 extreme corners of the (i,j) grid bound all cell centres (affine function of i,j
    ! optimised over a box always peaks at a vertex); pad by the point-to-point radius of a
    ! single hexagonal cell to safely enclose the whole lattice.
    c1 = self % cellCentre(1, 1)
    c2 = self % cellCentre(self % sizeN(1), 1)
    c3 = self % cellCentre(1, self % sizeN(2))
    c4 = self % cellCentre(self % sizeN(1), self % sizeN(2))

    halfwidth3 = ZERO
    halfwidth3(self % plane(1)) = (max(c1(1),c2(1),c3(1),c4(1)) - min(c1(1),c2(1),c3(1),c4(1))) * HALF &
                                  + hw * TWO_SQRT3
    halfwidth3(self % plane(2)) = (max(c1(2),c2(2),c3(2),c4(2)) - min(c1(2),c2(2),c3(2),c4(2))) * HALF &
                                  + hw * TWO_SQRT3
    halfwidth3(self % axis)     = self % sizeN(3) * self % axialPitch * HALF

    call tempDict % init(4)
    call tempDict % store('type', 'box')
    call tempDict % store('id', 1)
    call tempDict % store('origin', [ZERO, ZERO, ZERO])
    call tempDict % store('halfwidth', halfwidth3)
    call self % outline % init(tempDict)
    call tempDict % kill()

    ! Construct fill array
    call dict % get(tempI, 'map')

    ! Ensure size matches sizeN
    if (size(tempI) /= product(self % sizeN)) call fatalError(Here, &
            'Lattice map size not equal to size implied by shape. Respectively: '//&
            numToChar(size(tempI))//' '//numToChar(product(self % sizeN)))

    ! Flip array up-down for more natural input (as in latUniverse)
    ! Reshape into rank 2 array
    tempMap = reshape(tempI, [self % sizeN(1), self % sizeN(2) * self % sizeN(3)])
    N = size(tempMap, 2)
    do i = 1, N/2
      call swap(tempMap(:,i), tempMap(:,N - i + 1))
    end do

    ! Find background fill and change to tempMap to uniID
    tempMap = -tempMap
    call dict % get(name, 'padMat')
    outFill = charToFill(name, mats, Here)

    ! Build fill array
    allocate(fill(self % outLocalID))
    N = size(tempMap, 1)
    do j = 1, size(tempMap, 2)
      do i = 1, N
        fill(i + (j-1) * N) = tempMap(i, j)
      end do
    end do
    fill(self % outLocalID) = outFill
    deallocate(tempI, tempMap)

    ! Check whether there is an offset map
    if (dict % isPresent('offsetMap')) then

      if (.not. self % offset) call fatalError(Here, 'Cannot have both an offset map '//&
              'and no offset.')

      call dict % get(tempI, 'offsetMap')

      ! Ensure size matches sizeN
      if (size(tempI) /= product(self % sizeN)) call fatalError(Here, &
            'Offset map size not equal to size implied by shape. Respectively: '//&
            numToChar(size(tempI))//' '//numToChar(product(self % sizeN)))

      ! Flip array up-down for more natural input
      ! Reshape into rank 2 array
      tempMap = reshape(tempI, [self % sizeN(1), self % sizeN(2) * self % sizeN(3)])
      N = size(tempMap, 2)
      do i = 1, N/2
        call swap(tempMap(:,i), tempMap(:,N - i + 1))
      end do

      allocate(self % offsetMap(product(self % sizeN) + 1))
      N = size(tempMap, 1)
      do j = 1, size(tempMap, 2)
        do i = 1, N
          val = tempMap(i,j)
          self % offsetMap(i + (j-1) * N) = val

          ! Check that the entries are valid
          if (val /= local .and. val /= noOffset) call fatalError(Here,&
                  'Invalid entry to the offset map. Must be one of '//numToChar(local)//&
                  ' '//numToChar(noOffset)//'. Contains: '//numToChar(val))
        end do
      end do
      ! Add an entry for the padMat
      self % offsetMap(self % outLocalID) = noOffset

    end if

  end subroutine init

  !!
  !! Return the Cartesian (in-plane) position of the centre of cell (i,j)
  !!
  !! Args:
  !!   i [in] -> 1st in-plane oblique co-ordinate
  !!   j [in] -> 2nd in-plane oblique co-ordinate
  !!
  !! Result:
  !!   Cartesian position of the cell centre, in the local (plane) frame
  !!
  !! Note: valid (and used) for i,j outside the [1,N1] x [1,N2] range too, since the
  !!   underlying hexagonal tiling is infinite.
  !!
  pure function cellCentre(self, i, j) result(c)
    class(hexLatUniverse), intent(in) :: self
    integer(shortInt), intent(in)     :: i, j
    real(defReal), dimension(2)       :: c

    c = (i - (self % sizeN(1) + 1) * HALF) * self % a1 + &
        (j - (self % sizeN(2) + 1) * HALF) * self % a2

  end function cellCentre

  !!
  !! Find the (i,j) of the hexagonal cell (in the infinite tiling) whose centre is
  !! closest to the in-plane position rl
  !!
  !! Args:
  !!   rl [in] -> Position in the hexagonal plane (local frame)
  !!   i [out] -> 1st in-plane oblique co-ordinate of the closest cell
  !!   j [out] -> 2nd in-plane oblique co-ordinate of the closest cell
  !!
  pure subroutine nearestIJ(self, rl, i, j)
    class(hexLatUniverse), intent(in)        :: self
    real(defReal), dimension(2), intent(in)  :: rl
    integer(shortInt), intent(out)           :: i, j
    real(defReal), dimension(2)              :: latCoord, centre
    real(defReal)                            :: realI, realJ, dist2, minDist2
    integer(shortInt)                        :: baseI, baseJ, ti, tj

    ! Initial guess of the oblique (axial) co-ordinates
    latCoord = matmul(self % Ainv, rl)
    realI = latCoord(1) + (self % sizeN(1) + 1) * HALF
    realJ = latCoord(2) + (self % sizeN(2) + 1) * HALF
    baseI = nint(realI)
    baseJ = nint(realJ)

    ! Correct the guess: check the 3x3 neighbourhood for the truly closest hexagon centre
    ! (Simple rounding in oblique co-ordinates does not, in general, find the closest centre)
    minDist2 = INF
    i = baseI
    j = baseJ
    do ti = baseI - 1, baseI + 1
      do tj = baseJ - 1, baseJ + 1
        centre = self % cellCentre(ti, tj)
        dist2 = sum((rl - centre)**2)

        if (dist2 < minDist2) then
          minDist2 = dist2
          i = ti
          j = tj
        end if
      end do
    end do

  end subroutine nearestIJ

  !!
  !! Find local cell ID given a point
  !!
  !! See universe_inter for details.
  !!
  subroutine findCell(self, localID, cellIdx, r, u)
    class(hexLatUniverse), intent(inout)    :: self
    integer(shortInt), intent(out)          :: localID
    integer(shortInt), intent(out)          :: cellIdx
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    real(defReal), dimension(2)             :: rl, ul, centre, rl_local, n
    real(defReal)                           :: dist, proj, maxProj
    real(defReal)                           :: r_bar_ax
    integer(shortInt)                       :: i, j, k, f, bestFace, bestSign

    rl = r(self % plane)
    ul = u(self % plane)

    ! Find the (i,j) of the closest hexagonal cell in the infinite tiling
    call self % nearestIJ(rl, i, j)

    ! Identify which face (if any) is being crossed by the particle, including direction.
    ! This tie-break (as in hexagon_class % explicitBC) is necessary when the point sits on
    ! a corner between more than one cell.
    centre   = self % cellCentre(i, j)
    rl_local = rl - centre

    maxProj  = -INF
    bestFace = 0
    bestSign = 0
    do f = 1, 3
      n = self % faceNormal(f,:)
      dist = dot_product(rl_local, n)

      if (abs(dist) >= self % halfwidth - SURF_TOL) then
        proj = dot_product(ul, n * sign(ONE, dist))

        if (proj > maxProj) then
          maxProj  = proj
          bestFace = f
          bestSign = int(sign(ONE, dist))
        end if
      end if
    end do

    if (bestFace > 0 .and. maxProj > ZERO) then
      if (bestSign > 0) then
        i = i + self % faceOffset(bestFace, 1)
        j = j + self % faceOffset(bestFace, 2)
      else
        i = i - self % faceOffset(bestFace, 1)
        j = j - self % faceOffset(bestFace, 2)
      end if
    end if

    ! Find axial layer
    k = floor((r(self % axis) - self % axialCorner) / self % axialPitch) + 1
    r_bar_ax = r(self % axis) - self % axialCorner - (k - HALF) * self % axialPitch

    if (abs(r_bar_ax) > self % axialPitch * HALF - SURF_TOL .and. &
        r_bar_ax * u(self % axis) > ZERO) then
      if (u(self % axis) < ZERO) then
        k = k - 1
      else
        k = k + 1
      end if
    end if

    ! Set localID & cellIdx
    if (i < 1 .or. i > self % sizeN(1) .or. &
        j < 1 .or. j > self % sizeN(2) .or. &
        k < 1 .or. k > self % sizeN(3)) then ! Point is outside the lattice
      localID = self % outLocalID

    else
      localID = i + self % sizeN(1) * (j - 1 + self % sizeN(2) * (k - 1))

    end if
    cellIdx = 0

  end subroutine findCell

  !!
  !! Return distance to the next boundary between local cells in the universe
  !!
  !! See universe_inter for details.
  !!
  subroutine distance(self, d, surfIdx, coords)
    class(hexLatUniverse), intent(inout) :: self
    real(defReal), intent(out)           :: d
    integer(shortInt), intent(out)       :: surfIdx
    type(coord), intent(in)              :: coords
    real(defReal), dimension(2)          :: rl_local, ul, n
    real(defReal)                        :: uProj, rProj, a_far, test_d, r_bar_ax, centre_ax
    integer(shortInt), dimension(3)      :: ijk
    integer(shortInt)                    :: f, i, j, k

    ! Catch case if particle is outside the lattice
    if (coords % localID == self % outLocalID) then

      if (self % outline % evaluate(coords % r) > ZERO) then
        ! Genuinely outside the (padded) bounding box: safe to take one large jump to its
        ! surface, since the box is guaranteed to enclose the whole lattice.
        surfIdx = OUTLINE_SURF
        d = self % outline % distance(coords % r, coords % dir)
        return
      end if

      ! Otherwise: still background, but already inside the bounding box -- i.e. in the gap
      ! between the (padded) box and the true edge of the hexagonal tiling. Fall through to
      ! the standard per-cell face/axial test below, using the (i,j,k) of the tiling cell
      ! closest to the current position (even though it lies outside the valid map range),
      ! so the returned distance can never overshoot a genuine entry into the lattice.
      call self % nearestIJ(coords % r(self % plane), i, j)
      k = floor((coords % r(self % axis) - self % axialCorner) / self % axialPitch) + 1

    else
      ijk = get_ijk(coords % localID, self % sizeN)
      i = ijk(1); j = ijk(2); k = ijk(3)

    end if

    ! Position & direction wrt the centre of the current hexagonal cell
    rl_local = coords % r(self % plane) - self % cellCentre(i, j)
    ul       = coords % dir(self % plane)

    d = INF
    surfIdx = 0

    ! Loop over the 3 pairs of parallel hexagon edges
    do f = 1, 3
      n = self % faceNormal(f,:)
      uProj = dot_product(ul, n)
      rProj = dot_product(rl_local, n)

      if (uProj /= ZERO) then
        a_far  = sign(self % halfwidth, uProj)
        test_d = (a_far - rProj) / uProj
      else
        test_d = INF
      end if

      if (test_d < d) then
        d = test_d
        if (uProj > ZERO) then
          surfIdx = -(2 * f - 1)
        else
          surfIdx = -(2 * f)
        end if
      end if
    end do

    ! Axial direction
    centre_ax = self % axialCorner + (k - HALF) * self % axialPitch
    r_bar_ax  = coords % r(self % axis) - centre_ax
    uProj     = coords % dir(self % axis)

    if (uProj /= ZERO) then
      a_far  = sign(self % axialPitch * HALF, uProj)
      test_d = (a_far - r_bar_ax) / uProj
    else
      test_d = INF
    end if

    if (test_d < d) then
      d = test_d
      if (uProj > ZERO) then
        surfIdx = AX_MAX
      else
        surfIdx = AX_MIN
      end if
    end if

    ! Cap distance value
    d = max(ZERO, d)
    d = min(INF, d)

  end subroutine distance

  !!
  !! Cross between local cells
  !!
  !! See universe_inter for details.
  !!
  subroutine cross(self, coords, surfIdx)
    class(hexLatUniverse), intent(inout) :: self
    type(coord), intent(inout)           :: coords
    integer(shortInt), intent(in)        :: surfIdx

    call self % findCell(coords % localID, coords % cellIdx, coords % r, coords % dir)

  end subroutine cross

  !!
  !! Return offset for the current cell
  !!
  !! See universe_inter for details.
  !!
  function cellOffset(self, coords) result (offset)
    class(hexLatUniverse), intent(in) :: self
    type(coord), intent(in)           :: coords
    real(defReal), dimension(3)       :: offset
    logical(defBool)                  :: doOffset
    integer(shortInt), dimension(3)   :: ijk
    real(defReal), dimension(2)       :: centre

    if (allocated(self % offsetMap)) then
      doOffset = self % offsetMap(coords % localID) == local
    else
      doOffset = self % offset
    end if

    if (doOffset .and. coords % localID /= self % outLocalID) then
      ijk = get_ijk(coords % localID, self % sizeN)
      centre = self % cellCentre(ijk(1), ijk(2))

      offset = ZERO
      offset(self % plane) = centre
      offset(self % axis)  = self % axialCorner + (ijk(3) - HALF) * self % axialPitch

    else
      offset = ZERO

    end if

  end function cellOffset

  !!
  !! Return normal for the given surfIdx at a given point
  !!
  !! See universe_inter for details.
  !!
  function getNormal(self, surfIdx, coords) result (normal)
    class(hexLatUniverse), intent(in) :: self
    integer(shortInt), intent(in)     :: surfIdx
    type(coord), intent(in)           :: coords
    real(defReal), dimension(3)       :: normal
    integer(shortInt)                 :: f
    character(100), parameter :: Here = 'getNormal (hexLatUniverse_class.f90)'

    normal = ZERO

    select case(surfIdx)
      case(FACE1_POS, FACE1_NEG, FACE2_POS, FACE2_NEG, FACE3_POS, FACE3_NEG)
        f = (abs(surfIdx) + 1) / 2
        normal(self % plane) = self % faceNormal(f,:)
        if (mod(abs(surfIdx), 2) == 0) normal(self % plane) = -normal(self % plane)

      case(AX_MIN)
        normal(self % axis) = -ONE

      case(AX_MAX)
        normal(self % axis) = ONE

      case(OUTLINE_SURF)
        normal = self % outline % normal(coords % r, coords % dir)

      case default
        call fatalError(Here, 'Unrecognised surfIdx: '//numToChar(surfIdx))

    end select

  end function getNormal

  !!
  !! Return to uninitialised state
  !!
  subroutine kill(self)
    class(hexLatUniverse), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Local
    self % axis        = 0
    self % plane        = 0
    self % flatAxis     = 0
    self % pointAxis    = 0
    self % halfwidth    = ZERO
    self % axialPitch   = ZERO
    self % sizeN        = 0
    self % a1           = ZERO
    self % a2           = ZERO
    self % Ainv         = ZERO
    self % faceNormal   = ZERO
    self % faceOffset   = 0
    self % axialCorner  = ZERO
    call self % outline % kill()
    self % outLocalID   = 0
    self % offset       = .true.
    if (allocated(self % offsetMap)) deallocate(self % offsetMap)

  end subroutine kill

  !!
  !! Generate ijk from localID and shape
  !!
  !! Args:
  !!   localID [in] -> Local id of the cell between 1 and product(sizeN)
  !!   sizeN [in]   -> Number of cells in each direction: (N1, N2, N3)
  !!
  !! Result:
  !!   Array ijk which has integer position in each direction
  !!
  pure function get_ijk(localID, sizeN) result(ijk)
    integer(shortInt), intent(in)               :: localID
    integer(shortInt), dimension(3), intent(in) :: sizeN
    integer(shortInt), dimension(3)             :: ijk
    integer(shortInt)                           :: temp, base

    temp = localID - 1

    base = temp / sizeN(1)
    ijk(1) = temp - sizeN(1) * base + 1

    temp = base
    base = temp / sizeN(2)
    ijk(2) = temp - sizeN(2) * base + 1

    ijk(3) = base + 1

  end function get_ijk

end module hexLatUniverse_class
