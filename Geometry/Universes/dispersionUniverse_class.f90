module dispersionUniverse_class

  use numPrecision
  use universalVariables,     only : UNDEF_MAT, OVERLAP_MAT, INF
  use genericProcedures,      only : fatalError, numToChar
  use dictionary_class,       only : dictionary
  use coord_class,            only : coord
  use charMap_class,          only : charMap
  use surfaceShelf_class,     only : surfaceShelf
  use cellShelf_class,        only : cellShelf
  use sphere_class,           only : sphere
  use box_class,              only : box
  use cartesianLattice_class, only : cartesianLattice, X_MIN, X_MAX, Y_MIN, Y_MAX, &
                                          Z_MIN, Z_MAX
  use universe_inter,         only : universe, kill_super => kill, charToFill
  
  implicit none
  private

  ! Parameter defining the surface of the sphere bounding box
  integer(shortInt), parameter, public :: OUTLINE_SURF = Z_MAX -1

  ! Parameters defining the indices of various other fills.
  ! Background and outside are identical materials, but the IDs help cell searching
  ! and distance calculations: in the background, local spheres are checked. In the
  ! outside, no check is made for sphere intersections.
  integer(shortInt), parameter, public :: BACKGROUND_IDX = 1, OUTSIDE_IDX = 2,&
                                          OVERLAP_IDX = 3

  ! Parameter which oversizes the lattice by a small amount
  real(defReal), parameter, private :: OVER_SIZE = 1.0001_defReal

  !!
  !! Representation of a universe containing spherical inclusions (e.g., pebbles) in a matrix.
  !!
  !! The universe is created given a list of points, radii, and contents (materials or universes).
  !! A background material or universe is also required.
  !!
  !! The input consists of a file, containing the radii, origin, and fill information. This
  !! is very similar to Serpent's HTGR universe. The syntax for this is:
  !!   radius, x, y, z, localFill
  !! For example:
  !!   7.0 -1 5.3 2.222 pyC
  !!   4 14.2 46 -10.0 u<3>
  !!
  !! Each local cell in the universe corresponds to either the background material or a
  !! spherical inclusion.
  !! Two extra local cells are always defined inside the dispersionUniverse. One has OVERLAP_MAT
  !! (overlapping cells) filling. The other has the background fill, but is OUTSIDE the dispersion 
  !! region.
  !! If a position is in more than one user-defined cell, it is in the overlapping cell. 
  !! This exists to enable plotting of geometry without fatalErrors.
  !! However, overlapping cells can only be detected if a less optimal cell search is enabled
  !! using the checkOverlap flag. This is encouraged for plotting and debugging, but discouraged
  !! during transport.
  !! The OUTSIDE region exists to avoid searching for sphere intersections in case the particle is
  !! in this region. This is distinct from the BACKGROUND region where spheres may be nearby.
  !!
  !! Sample Input Dictionary:
  !!   uni { type dispersionUniverse;
  !!         id 7;
  !!         file /path/to/file.txt;
  !!         background graphite;
  !!         # rotation (23.0 0.0 0.0); #
  !!         # checkOverlap 0;          #
  !!       }
  !!
  !! Private Members:
  !!   sphere           -> array of sphere surfaces which define the dispersion
  !!   outline          -> box bounding all spheres
  !!   accelStruct      -> lattice which contains all spheres and provides rapid 
  !!                       cell finding and local distance calculations
  !!   sphereInGridList -> array of sphere indices which are in a given lattice grid cell.
  !!                       Used with offsetToGridCell to compactly store sphere/grid info.
  !!   offsetToGridCell -> array of offsets between grid cell entries in the sphereInGridList
  !!   avgRadius        -> average radius of all spheres. Used to estimate the grid pitch.
  !!   checkOverlap     -> flag determining whether to perform more expensive but robust
  !!                       sphere searches.
  !!   gridSize         -> Number of elements in the search grid
  !!   nSpheres         -> Number of spheres in the universe
  !!
  !! Interface:
  !!   universe interface
  !!
  type, public, extends(universe) :: dispersionUniverse
    private
    type(sphere), dimension(:), allocatable      :: spheres
    type(box)                                    :: outline
    type(cartesianLattice)                       :: accelStruct
    integer(shortInt), dimension(:), allocatable :: sphereInGridList
    integer(shortInt), dimension(:), allocatable :: offsetToGridCell
    real(defReal)                                :: avgRadius = ZERO
    logical(defBool)                             :: checkOverlap = .false.
    integer(shortInt)                            :: gridSize = 0
    integer(shortInt)                            :: nSpheres = 0
  contains
    ! Superclass procedures
    procedure :: init
    procedure :: kill
    procedure :: findCell
    procedure :: distance
    procedure :: cross
    procedure :: cellOffset
    procedure :: getNormal
    ! Local procedures
    procedure, private :: readFile
    procedure, private :: buildAccelerationStructure
  end type dispersionUniverse

contains

  !!
  !! Initialise Universe
  !!
  !! See universe_inter for details.
  !!
  subroutine init(self, fill, dict, cells, surfs, mats)
    class(dispersionUniverse), intent(inout)                  :: self
    integer(shortInt), dimension(:), allocatable, intent(out) :: fill
    class(dictionary), intent(in)                             :: dict
    type(cellShelf), intent(inout)                            :: cells
    type(surfaceShelf), intent(inout)                         :: surfs
    type(charMap), intent(in)                                 :: mats
    character(nameLen)                                        :: background
    integer(shortInt)                                         :: backgroundFill
    character(pathLen)                                        :: fileName
    character(100), parameter :: Here = 'init (dispersionUniverse_class.f90)'

    ! Setup the base class
    ! With: id, origin rotations...
    call self % setupBase(dict)

    ! Read the background material
    call dict % get(background, 'background')
    backgroundFill = charToFill(background, mats, Here)

    ! Get the file
    call dict % get(fileName, 'file')

    ! Open the file and check the number and validity of entries
    call self % readFile(fileName, mats, backgroundFill, fill)

    ! Build the acceleration structure
    call self % buildAccelerationStructure()

    ! Check for overlaps?
    call dict % getOrDefault(self % checkOverlap, 'checkOverlap', .false.)

  end subroutine init

  !!
  !! Read file containing positions, radii, and fill information
  !! This should be in the format:
  !!   radius x y z fill
  !!
  !! FatalError if file does not exist or is invalid
  !!
  subroutine readFile(self, fileName, mats, backgroundFill, fill)
    class(dispersionUniverse), intent(inout)                  :: self
    character(pathLen), intent(in)                            :: fileName
    type(charMap), intent(in)                                 :: mats
    integer(shortInt), intent(in)                             :: backgroundFill
    integer(shortInt), dimension(:), allocatable, intent(out) :: fill
    integer(shortInt)                                         :: i, N
    real(defReal), dimension(3)                               :: rMin, rMax
    real(defReal)                                             :: radius, x, y, z
    character(nameLen)                                        :: localFill
    type(dictionary)                                          :: tempDict
    logical(defBool)                                          :: exists
    character(100), parameter :: Here = 'readFile (dispersionUniverse_class.f90)'

    ! Check the file exists
    inquire(file=fileName, exist=exists)
    if (.not. exists) then
      call fatalError(Here, 'File does not exist: '//fileName)
    end if

    ! Open the file and scan each line
    open(unit=10, file=fileName, status='old', action='read', iostat=i)
    if (i /= 0) call fatalError(Here, 'Could not open file: '//fileName)

    ! Read each line and check the maximum values
    N = 0
    rMin = [INF, INF, INF]
    rMax = [-INF, -INF, -INF]
    do
      read(10, *, iostat=i) radius, x, y, z, localFill
      if (i /= 0) exit
      N = N + 1

      if (radius <= ZERO) then
        call fatalError(Here, 'Radius must be positive in file: '//fileName)
      end if

      rMin = min(rMin, [x - radius, y - radius, z - radius])
      rMax = max(rMax, [x + radius, y + radius, z + radius])

    end do
    close(10)

    if (N == 0) call fatalError(Here,'Universe definition contains no inclusions.')
    self % nSpheres = N

    ! Build outline box
    ! Slightly oversize it to avoid the lattice missing points on the edge
    call tempDict % init(4)
    call tempDict % store('type', 'box')
    call tempDict % store('id', 1)
    call tempDict % store('origin', (rMax + rMin) * HALF)
    call tempDict % store('halfwidth', abs(rMax-rMin) * HALF * OVER_SIZE)
    call self % outline % init(tempDict)
    call tempDict % kill()

    if (allocated(self % spheres)) deallocate(self % spheres)
    allocate(self % spheres(self % nSpheres))
    allocate(fill(self % nSpheres + OVERLAP_IDX))  ! +3 for background, outside, OVERLAP cells

    ! Do the second file pass
    open(unit=10, file=fileName, status='old', action='read', iostat=i)
    do i = 1, self % nSpheres
      read(10, *) radius, x, y, z, localFill
      call self % spheres(i) % build(1, [x, y, z], radius)
      fill(i) = charToFill(localFill, mats, Here)
      self % avgRadius = self % avgRadius + radius
    end do
    close(10)
    fill(self % nSpheres + BACKGROUND_IDX) = backgroundFill
    fill(self % nSpheres + OUTSIDE_IDX) = backgroundFill
    fill(self % nSpheres + OVERLAP_IDX) = OVERLAP_MAT
    self % avgRadius = self % avgRadius / real(self % nSpheres, defReal)

  end subroutine readFile

  !!
  !! Overlays a regular grid on the universe to accelerate cell searches. The grid is uniform in
  !! each dimension. A search over the spheres is performed to determine which overlap with a given
  !! grid cell.
  !!
  subroutine buildAccelerationStructure(self)
    class(dispersionUniverse), intent(inout)     :: self
    integer(shortInt), dimension(3)              :: n
    real(defReal)                                :: pitch, radius
    type(dictionary)                             :: tempDict
    integer(shortInt)                            :: s, i, j, k, idx, pos
    integer(shortInt), dimension(3)              :: ijkMin, ijkMax
    integer(shortInt), dimension(:), allocatable :: counts, fillPos
    real(defReal), dimension(6)                  :: box, outside
    real(defReal), dimension(3)                  :: origin, pitch3
    character(100), parameter :: Here = 'buildAccelerationStructure (dispersionUniverse_class.f90)'
  
    if (self % avgRadius <= ZERO) call fatalError(Here,&
        "Invalid average sphere radius: "//numToChar(self % avgRadius))

    ! Create the cartesian lattice
    pitch = self % avgRadius * TWO
    outside = self % outline % boundingBox()
    n = ceiling([outside(4:6) - outside(1:3)] / pitch)

    ! Recompute the pitch
    pitch3 = (outside(4:6) - outside(1:3)) / n
    call tempDict % init(2)
    call tempDict % store('pitch', pitch3)
    call tempDict % store('shape', n)
    call tempDict % store('origin', HALF * (outside(1:3)+outside(4:6)))
    call self % accelStruct % init(tempDict)

    self % gridSize = product(self % accelStruct % getSize())

    ! Do a first pass through each to determine how many intersections spheres have with 
    ! a given grid element.
    allocate(counts(self % gridSize))
    counts = 0

    do s = 1, self % nSpheres

      ! Check bounding box intersections with the acceleration structure
      box = self % spheres(s) % boundingBox()
      origin = self % spheres(s) % getOrigin()
      radius = self % spheres(s) % getRadius()

      ijkMin = self % accelStruct % findIJK(box(1:3))
      ijkMax = self % accelStruct % findIJK(box(4:6))

      ! Catch the indices being out of bounds: the outline should have been sized
      ! to prevent this!
      if (any(ijkMin == 0) .or. any(ijkMax == 0)) then
        call fatalError(Here,'Inclusion with radius '//numToChar(radius)//' and origin '//&
            numToChar(origin)//' has a box outside of the outline, which is: '//&
            numToChar(outside))
      end if
      
      ! Check all cells from ijkMin to ijkMax for sphere intersection
      do k = ijkMin(3), ijkMax(3)
        do j = ijkMin(2), ijkMax(2)
          do i = ijkMin(1), ijkMax(1)
            idx = self % accelStruct % findCell([i, j, k])
            if (sphereOverlap(origin, radius, self % accelStruct % getBox(idx))) then
              counts(idx) = counts(idx) + 1
            end if
          end do
        end do
      end do

    end do

    ! Information is stored in a CSR format:
    ! In order of grid cells, the spheres which intersect that cell are recorded in an array.
    ! Subsequently, the offset of the first sphere for each grid cell is recorded in a second array.
    if (any(counts > 0)) then
      if (allocated(self % sphereInGridList)) deallocate(self % sphereInGridList)
      if (allocated(self % offsetToGridCell)) deallocate(self % offsetToGridCell)
      allocate(self % sphereInGridList(sum(counts)))
      allocate(self % offsetToGridCell(self % gridSize + 1))

      self % sphereInGridList = 0
      self % offsetToGridCell = 0

      self % offsetToGridCell(1) = 1
      
      if (size(self % offsetToGridCell) > 1) then
        do i = 1, self % gridSize
          self % offsetToGridCell(i+1) = self % offsetToGridCell(i) + counts(i)
        end do
      end if

      ! fillPos tracks the next free write slot for each cell, starts at offsets
      allocate(fillPos(self % gridSize))
      fillPos = self % offsetToGridCell(1:self % gridSize)

      do s = 1, self % nSpheres

        ! Check bounding box intersections with the acceleration structure
        box = self % spheres(s) % boundingBox()
        origin = self % spheres(s) % getOrigin()
        radius = self % spheres(s) % getRadius()

        ijkMin = self % accelStruct % findIJK(box(1:3))
        ijkMax = self % accelStruct % findIJK(box(4:6))
      
        ! Check all cells from ijkMin to ijkMax for sphere intersection
        do k = ijkMin(3), ijkMax(3)
          do j = ijkMin(2), ijkMax(2)
            do i = ijkMin(1), ijkMax(1)
              idx = self % accelStruct % findCell([i, j, k])
              if (sphereOverlap(origin, radius, self % accelStruct % getBox(idx))) then
                pos = fillPos(idx)
                self % sphereInGridList(pos) = s
                fillPos(idx) = pos + 1
              end if
            end do
          end do
        end do

      end do

      deallocate(fillPos)

    else

      ! Allocate all offsets to 1
      ! Should avoid looping over spheres in findCell due to
      ! the resulting do being from 1 to 0, which doesn't execute
      allocate(self % offsetToGridCell(self % gridSize + 1))
      self % offsetToGridCell = 1
      allocate(self % sphereInGridList(0))

    end if

  end subroutine buildAccelerationStructure

  !!
  !! Find local cell ID given a point
  !!
  !! See universe_inter for details.
  !!
  subroutine findCell(self, localID, cellIdx, r, u)
    class(dispersionUniverse), intent(inout) :: self
    integer(shortInt), intent(out)           :: localID
    integer(shortInt), intent(out)           :: cellIdx
    real(defReal), dimension(3), intent(in)  :: r
    real(defReal), dimension(3), intent(in)  :: u
    integer(shortInt)                        :: found, i, s, idx

    ! cellIdx is always zero - no cellShelf cells in this universe
    cellIdx = 0

    ! Count how many valid spheres have been found. Guard against overlaps
    found = 0

    ! Find index in the acceleration structure
    idx = self % accelStruct % findCell(r, u)

    ! If the point is outside the acceleration structure, it is in the background cell
    if (idx > self % gridSize) then
      localID = size(self % spheres) + OUTSIDE_IDX
      return
    end if

    ! Otherwise, need to search any constituent spheres which overlap with the grid cell
    do i = self % offsetToGridCell(idx), self % offsetToGridCell(idx + 1) - 1
      s = self % sphereInGridList(i)
      if (.not. self % spheres(s) % halfspace(r, u)) then
        localID = s
        if (self % checkOverlap) then
          found = found + 1
        else
          return
        end if
      end if
    end do

    if (self % checkOverlap) then
      if (found == 1) then
        return

      ! Return OVERLAP_MAT if more than one sphere is valid
      else if (found > 1) then
        localID = size(self % spheres) + OVERLAP_IDX
        return
      end if
    end if
      
    ! If not in any of the spheres, must be in the background
    localID = self % nSpheres + BACKGROUND_IDX

  end subroutine findCell

  !!
  !! Return distance to the next boundary between local cells in the universe
  !!
  !! See universe_inter for details.
  !!
  !! Errors:
  !!   fatalError if in OVERLAP cell
  !!
  subroutine distance(self, d, surfIdx, coords)
    class(dispersionUniverse), intent(inout) :: self
    real(defReal), intent(out)               :: d
    integer(shortInt), intent(out)           :: surfIdx
    type(coord), intent(in)                  :: coords
    integer(shortInt)                        :: localID, idx, i, s
    real(defReal)                            :: dTest
    real(defReal), dimension(3)              :: r, u
    character(100), parameter :: Here = 'distance (dispersionUniverse_class.f90)'

    localID = coords % localID

    r = coords % r
    u = coords % dir

    ! Catch case if particle is outside the lattice
    if (coords % localID == self % nSpheres + OUTSIDE_IDX) then
      surfIdx = OUTLINE_SURF
      d = self % outline % distance(r, u)
      return

    ! Also catch if particle is in overlapping material
    else if (coords % localID == self % nSpheres + OVERLAP_IDX) then
      call fatalError(Here, 'Particle is in an overlapping local cell. Position: '&
          //numToChar(coords % r))

    end if

    ! Otherwise compute distance to grid and sphere boundaries

    ! Find index in the acceleration structure
    idx = self % accelStruct % findCell(r, u)
    call self % accelStruct % distance(d, surfIdx, idx, r, u)

    ! Test distance to each local sphere
    do i = self % offsetToGridCell(idx), self % offsetToGridCell(idx + 1) - 1
      
      s = self % sphereInGridList(i)
      dTest = self % spheres(s) % distance(r, u)
      
      ! The surfIdx -s + OUTLINE_SURF is because 1) negative surfIdxs show that surfaces are not
      ! on the surfShell but instead part of a universe and 2) surfIdx -7 is the outline
      ! so -s+OUTLINE_SURF is the index of sphere s
      if (dTest < d) then
        d = dTest
        surfIdx = -s + OUTLINE_SURF
      end if
    
    end do

  end subroutine distance

  !!
  !! Cross between local cells
  !!
  !! See universe_inter for details.
  !!
  subroutine cross(self, coords, surfIdx)
    class(dispersionUniverse), intent(inout) :: self
    type(coord), intent(inout)               :: coords
    integer(shortInt), intent(in)            :: surfIdx

    call self % findCell(coords % localID, coords % cellIdx, coords % r, coords % dir)

  end subroutine cross

  !!
  !! Return offset for the current cell
  !!
  !! See universe_inter for details.
  !!
  function cellOffset(self, coords) result (offset)
    class(dispersionUniverse), intent(in) :: self
    type(coord), intent(in)               :: coords
    real(defReal), dimension(3)           :: offset

    if (coords % localID <= self % nSpheres) then
      offset = self % spheres(coords % localID) % getOrigin()
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
    class(dispersionUniverse), intent(in) :: self
    integer(shortInt), intent(in)         :: surfIdx
    type(coord), intent(in)               :: coords
    real(defReal), dimension(3)           :: normal
    integer(shortInt)                     :: idx
    character(100), parameter :: Here = 'getNormal (dispersionUniverse_class.f90)'
    
    ! Internal boundary within the acceleration structure
    if (surfIdx > OUTLINE_SURF .and. surfIdx <= X_MIN) then
      normal = self % accelStruct % getNormal(surfIdx)

    ! Hits the outline surface
    else if (surfIdx == OUTLINE_SURF) then
      normal = self % outline % normal(coords % r, coords % dir)

    ! Hits one of the spheres
    else if (surfIdx < OUTLINE_SURF .and. surfIdx >= (OUTLINE_SURF - self % nSpheres)) then
      idx = abs(surfIdx - OUTLINE_SURF)
      normal = self % spheres(idx) % normal(coords % r, coords % dir)
    
    else
      call fatalError(Here, "Invalid surfIdx: "//numToChar(surfIdx))
    end if

  end function getNormal

  !!
  !! Return to uninitialised state
  !!
  subroutine kill(self)
    class(dispersionUniverse), intent(inout) :: self
    integer(shortInt)                        :: s

    ! Superclass
    call kill_super(self)

    ! Local
    if(allocated(self % spheres)) then
      do s = 1, self % nSpheres
        call self % spheres(s) % kill()
      end do
      deallocate(self % spheres)
    end if
    call self % outline % kill()
    call self % accelStruct % kill()
    if(allocated(self % sphereInGridList)) deallocate(self % sphereInGridList)
    if(allocated(self % offsetToGridCell)) deallocate(self % offsetToGridCell)
    self % avgRadius = ZERO
    self % gridSize = 0
    self % nSpheres = 0

  end subroutine kill

  pure function sphereOverlap(origin, radius, box) result(overlap)
    real(defReal), dimension(3), intent(in) :: origin
    real(defReal), intent(in)               :: radius
    real(defReal), dimension(6), intent(in) :: box
    logical(defBool)                        :: overlap
    real(defReal), dimension(3)             :: closestPoint

    ! Find the closest point in the box to the sphere's origin
    closestPoint = max(box(1:3), min(origin, box(4:6)))

    ! Check if the distance from the closest point to the sphere's origin is less than the radius
    overlap = sum((closestPoint - origin)**2) <= radius * radius

  end function sphereOverlap

end module dispersionUniverse_class
