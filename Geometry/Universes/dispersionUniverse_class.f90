module dispersionUniverse_class

  use numPrecision
  use universalVariables, only : UNDEF_MAT, OVERLAP_MAT, NUDGE
  use genericProcedures,  only : fatalError, numToChar
  use dictionary_class,   only : dictionary
  use coord_class,        only : coord
  use charMap_class,      only : charMap
  use surfaceShelf_class, only : surfaceShelf
  use universe_inter,     only : universe, kill_super => kill
  
  implicit none
  private

  ! Parameter which identifies whether a cell in the acceleration structure contains
  ! no particles.
  integer(shortInt), parameter, public :: NO_PARTICLE = -1

  !!
  !! Representation of a universe containing spherical inclusions (e.g., pebbles) in a matrix.
  !!
  !! The universe is created given a list of points, radii, and contents (materials or universes).
  !! A background material or universe is also required.
  !!
  !! The input consists of a file, containing the position, radii, and fill information. This
  !! is very similar to Serpent's HTGR universe.
  !!
  !! Each local cell in the universe corresponds to a cell given by an ID
  !! Two extra local cells are always defined inside the dispersionUniverse with UNDEF_MAT 
  !! (undefined material) and OVERLAP_MAT (overlapping cells) filling.
  !! If position is not in any user-defined cell, it is in the undefined cell. If position
  !! is in more than one user-defined cell, it is in the overlapping cell. These both exist
  !! to enable plotting of geometry without fatalErrors.
  !! However, overlapping cells can only be detected if a less optimal cell search is enabled
  !! using the checkOverlap flag. This is encouraged for plotting and debugging, but discouraged
  !! during transport.
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
  !! Note:
  !!   - Cell overlaps are forbidden, but there is no check to find overlaps unless specified in input.
  !!
  !! Public Members:
  !!
  !!
  !! Interface:
  !!   universe interface
  !!
  type, public, extends(universe) :: dispersionUniverse
    type(sphere), dimension(:), allocatable    :: spheres
    type(box)                                  :: outline
    type(cartesianLattice)                     :: accelStruct
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
    procedure :: readFile
    procedure :: buildAccelerationStructure
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
    backgroundFill = charToFill(background)

    ! Get the file
    call dict % get(fileName, 'file')

    ! Open the file and check the number and validity of entries
    call self % readFile(fileName, fill, backgroundFill)

    ! Build the acceleration structure
    call self % buildAccelerationStructure()

    ! Check for overlaps?
    call dict % getOrDefault(self % checkOverlap, 'checkOverlap', .false.)

  end subroutine init

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
    integer(shortInt)                        :: i, tIdx
    integer(shortInt), dimension(:), allocatable :: array

    ! If being careful, only do exhaustive searches to check for overlap
    if (self % checkOverlap) then
      call self % findCellOverlap(localID, cellIdx, r, u)
      return
    end if

  end subroutine findCell

  !!
  !! Return distance to the next boundary between local cells in the universe
  !!
  !! See universe_inter for details.
  !!
  !! Errors:
  !!   fatalError if in UNDEFINED cell or OVERLAP cell
  !!
  subroutine distance(self, d, surfIdx, coords)
    class(dispersionUniverse), intent(inout) :: self
    real(defReal), intent(out)         :: d
    integer(shortInt), intent(out)     :: surfIdx
    type(coord), intent(in)            :: coords
    integer(shortInt)                  :: localID
    character(100), parameter :: Here = 'distance (dispersionUniverse_class.f90)'

    localID = coords % localID

    if (localID == size(self % cells) + 1) then
      call fatalError(Here, 'Particle is in undefined local cell. Local ID: '//numToChar(localID))
    elseif (localID == size(self % cells) + 2) then
      call fatalError(Here, 'Particle is in an overlapping local cell. Local ID: '//numToChar(localID))
    end if

    ! Calculate distance
    call self % cells(localID) % ptr % distance(d, surfIdx, coords % r, coords % dir)

  end subroutine distance

  !!
  !! Cross between local cells
  !!
  !! See universe_inter for details.
  !!
  !! Note: Introduces extra movement to the particle to push it over boundary
  !!   for more efficient search. Distance is NUGDE.
  !!
  subroutine cross(self, coords, surfIdx)
    class(dispersionUniverse), intent(inout) :: self
    type(coord), intent(inout)         :: coords
    integer(shortInt), intent(in)      :: surfIdx
    integer(shortInt)                  :: local0
    logical(defBool)                   :: foundNeighb

    ! Keep initial cell ID to add to neighbour lists
    local0 = coords % localID

    ! NUDGE position slightly forward to escape surface tolerance
    ! and avoid calculating normal and extra dot-products
    coords % r = coords % r + coords % dir * NUDGE

    ! Find cell
    ! First perform a neighbour search
    if (.not. self % checkOverlap) then
      call self % findCellNeighb(coords % localID, &
                           coords % cellIdx, &
                           coords % r,       &
                           coords % dir,     &
                           local0,           &
                           foundNeighb)

      if (foundNeighb) return
    end if
    
    ! If that failed, perform an exhaustive search and add the
    ! new cell to the neighbour list
    call self % findCell(coords % localID, &
                         coords % cellIdx, &
                         coords % r,       &
                         coords % dir)

    ! Ensure the cells can be added to the neighbour list
    if (coords % cellIdx < 1) return

    ! Add each cell to the other's neighbour list
    if (.not. self % checkOverlap) then
      call self % cells(local0) % neighb % add(coords % localID)
      call self % cells(coords % localID) % neighb % add(local0)
    end if

  end subroutine cross

  !!
  !! Return offset for the current cell
  !!
  !! See universe_inter for details.
  !!
  function cellOffset(self, coords) result (offset)
    class(dispersionUniverse), intent(in) :: self
    type(coord), intent(in)         :: coords
    real(defReal), dimension(3)     :: offset

    ! There is no cell offset
    offset = ZERO

  end function cellOffset
    
  !!
  !! Return normal for the given surfIdx at a given point
  !!
  !! See universe_inter for details.
  !!
  function getNormal(self, surfIdx, coords) result (normal)
    class(dispersionUniverse), intent(in) :: self
    integer(shortInt), intent(in)   :: surfIdx
    type(coord), intent(in)         :: coords
    real(defReal), dimension(3)     :: normal
    integer(shortInt)               :: cIdx
    character(100), parameter :: Here = 'getNormal (dispersionUniverse_class.f90)'

    ! Ensure that the current cell of coords has surfIdx as a component
    cIdx = coords % localID

    ! Ensure that there are this many cells
    if (cIdx > size(self % cells) .or. cIdx < 1) call fatalError(Here, &
           'cIdx is invalid: '//numToChar(cIdx))

    ! Ensure that the surfIdx is valid
    if (surfIdx < 1) call fatalError(Here, &
           'surfIdx is invalid: '//numToChar(surfIdx)) 

    normal = self % cells(cIdx) % ptr % getNormal(surfIdx, coords % r, coords % dir)   

  end function getNormal

  !!
  !! Return to uninitialised state
  !!
  subroutine kill(self)
    class(dispersionUniverse), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Local
    if(allocated(self % radii)) deallocate(self % radii)
    if(allocated(self % positions)) deallocate(self % positions)
    self % N = 0
    call self % outline % kill()

  end subroutine kill

  !!
  !! Read file containing positions, radii, and fill information
  !! This should be in the format:
  !!   radius x y z fill
  !!
  !! FatalError if file does not exist or is invalid
  !!
  subroutine readFile(self, fileName, backgroundFill, fill)
    class(dispersionUniverse), intent(inout)                  :: self
    character(pathLen), intent(in)                            :: fileName
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
    if (.not. fileExists(fileName)) then
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
      read(10, *, iostat=i) radius, x, y, z, fill
      if (i /= 0) exit
      N = N + 1

      if (radius <= 0.0) then
        call fatalError(Here, 'Radius must be positive in file: '//fileName)
      end if

      rMin = min(rMin, [x - radius, y - radius, z - radius])
      rMax = max(rMax, [x + radius, y + radius, z + radius])

    end do
    close(10)

    ! Build outline box
    call tempDict % init(4)
    call tempDict % store('type', 'box')
    call tempDict % store('id', 1)
    call tempDict % store('origin', (rMax + rMin) * HALF)
    call tempDict % store('halfwidth', abs(rMax-rMin) * HALF)
    call self % outline % init(tempDict)
    call tempDict % kill()

    allocate(self % spheres(N))
    allocate(fill(N+3))  ! +3 for background, UNDEF, and OVERLAP cells

    ! Do the second file pass
    open(unit=10, file=fileName, status='old', action='read', iostat=i)
    do i = 1, N
      read(10, *) radius, x, y, z, localFill
      call tempDict % init(4)
      call tempDict % store('type', 'sphere')
      call tempDict % store('id', i+1)
      call tempDict % store('origin', [x, y, z])
      call tempDict % store('radius', radius)
      call self % spheres(i) % init(tempDict)
      call tempDict % kill()
      fill(i) = charToFill(localFill)
      self % avgRadius = self % avgRadius + radius
    end do
    close(10)
    fill(N+1) = backgroundFill
    fill(N+2) = UNDEF_MAT
    fill(N+3) = OVERLAP_MAT
    self % avgRadius = self % avgRadius / real(N, defReal)

  end subroutine readFile

  !!
  !! Overlays a regular grid on the universe to accelerate cell searches. The grid is uniform in
  !! each dimension. A search over the spheres is performed to determine which overlap with a given
  !! grid cell.
  !!
  subroutine buildAccelerationStructure(self)
    class(dispersionUniverse), intent(inout) :: self
    integer(shortInt), dimension(3)          :: n
    real(defReal)                            :: pitch
    real(defReal), dimension(6)              :: box
    character(100), parameter :: Here = 'buildAccelerationStructure (dispersionUniverse_class.f90)'

    ! Create the cartesian lattice
    pitch = self % avgRadius * TWO
    box = self % outline % boundingBox()
    n = ceiling([box(4:6) - box(1:3)] / pitch)
    call tempDict % init(2)
    call tempDict % store('pitch', [pitch, pitch, pitch])
    call tempDict % store('shape', n)
    call self % accelStruct % init(tempDict)

    ! Do a first pass through each to determine how many intersections they have with 
    ! a given mesh element.
    allocate(counts(product(self % accelStruct % getSize())))
    do i = 1, size(self % spheres)

      ! Check bounding box intersections with the acceleration structure
      box = self % spheres(i) % boundingBox()
      ijkMin = self % accelStruct % getIJK(box(1:3))
      ijkMax = self % accelStruct % getIJK(box(4:6))
      call self % accelStruct % countIntersections(self % spheres(i))
    end do

    ! For each mesh element, allocate the correct number of intersections
    do i = 1, product(self % accelStruct % getSize())
      allocate(self % )
      call self % accelStruct % allocateIntersections(i)
    end do



    

   

  end subroutine buildAccelerationStructure

end module dispersionUniverse_class
