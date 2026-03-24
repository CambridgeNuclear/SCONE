module cellUniverse_class

  use numPrecision
  use universalVariables, only : UNDEF_MAT, OVERLAP_MAT, NUDGE
  use genericProcedures,  only : fatalError, numToChar
  use dictionary_class,   only : dictionary
  use coord_class,        only : coord
  use charMap_class,      only : charMap
  use linkedList_class,   only : linkedIntList
  use surfaceShelf_class, only : surfaceShelf
  use cell_inter,         only : cell
  use cellShelf_class,    only : cellShelf
  use universe_inter,     only : universe, kill_super => kill
  use omp_lib,            only : omp_test_lock, omp_unset_lock, &
                                 omp_init_lock, omp_destroy_lock, &
                                 omp_lock_kind
  implicit none
  private
  
  ! Determines how frequently to reorder the cell list
  integer(longInt), public, parameter :: reorderFreq = 10_longInt**6
  
  !!
  !! Local helper class to group cell data
  !!
  !! Public Members:
  !!   idx    -> cellIdx of the cell in cellShelf
  !!   ptr    -> Pointer to the cell
  !!   visits -> Number of times the cell has been found
  !!   neighb -> Neighbour list of adjacent cells
  !!
  type, private :: localCell
    integer(shortInt)    :: idx = 0
    class(cell), pointer :: ptr => null()
    integer(longInt)     :: visits = 0
    type(linkedIntList)  :: neighb
  end type localCell

  !!
  !! Representation of a universe via cells
  !!
  !! Each local cell in the universe corresponds to a cell given by an ID.
  !! Two extra local cells are always defined inside the cellUniverse with UNDEF_MAT 
  !! (undefined material) and OVERLAP_MAT (overlapping cells) filling.
  !! If position is not in any user-defined cell, it is in the undefined cell. If position
  !! is in more than one user-defined cell, it is in the overlapping cell. These both exist
  !! to enable plotting of geometry without fatalErrors.
  !! However, overlapping cells can only be detected if a less optimal cell search is enabled
  !! using the checkOverlap flag. This is encouraged for plotting and debugging, but discouraged
  !! during transport.
  !!
  !! The universe periodically reorders the cells by frequency of visits to accelerate cell
  !! search operations.
  !! Additionally, neighbour lists are built when the cross operation is invoked. These also
  !! accelerate cell searches.
  !!
  !! Sample Input Dictionary:
  !!   uni { type cellUniverse;
  !!         id 7;
  !!         # origin (2.0 0.0 0.0);    #
  !!         # rotation (23.0 0.0 0.0); #
  !!         # checkOverlap 0;          #
  !!         cells ( 1 3 4);         }
  !!
  !! Note:
  !!   - Local IDs are assigned in order as in definition. In the example above local id would map
  !!   to following cell  [localID: cellID] {1: 1, 2: 3, 3: 4, 4: UNDEF, 5: OVERLAP }
  !!   - Cell overlaps are forbidden, but there is no check to find overlaps unless specified in input.
  !!
  !! Public Members:
  !!   cells        -> Structure that stores cellIdx and pointers to the cells
  !!   searchIdx    -> Array of cell indices to search in order. Used to accelerate search.
  !!   reordering   -> Flag to allow searchIdx to be safely reordered in parallel
  !!   checkOverlap -> Flag to induce slower, thorough cell searches to detect overlaps
  !!
  !! Interface:
  !!   universe interface
  !!
  type, public, extends(universe) :: cellUniverse
    type(localCell), dimension(:), allocatable   :: cells
    integer(shortInt), dimension(:), allocatable :: searchIdx
    logical(defBool) :: reordering = .false.
    logical(defBool) :: checkOverlap = .false.
  contains
    ! Superclass procedures
    procedure :: init
    procedure :: kill
    procedure :: findCell
    procedure :: distance
    procedure :: cross
    procedure :: cellOffset
    
    ! Local procedures
    procedure :: findCellOverlap
    procedure :: findCellNeighb
    procedure :: reorderCells
  end type cellUniverse

contains

  !!
  !! Initialise Universe
  !!
  !! See universe_inter for details.
  !!
  subroutine init(self, fill, dict, cells, surfs, mats)
    class(cellUniverse), intent(inout)                        :: self
    integer(shortInt), dimension(:), allocatable, intent(out) :: fill
    class(dictionary), intent(in)                             :: dict
    type(cellShelf), intent(inout)                            :: cells
    type(surfaceShelf), intent(inout)                         :: surfs
    type(charMap), intent(in)                                 :: mats
    integer(shortInt), dimension(:), allocatable  :: cellTemp
    integer(shortInt)                             :: N, i
    character(100), parameter :: Here = 'init (cellUniverse_class.f90)'

    ! Setup the base class
    ! With: id, origin rotations...
    call self % setupBase(dict)

    ! Load Cells by ID
    call dict % get(cellTemp, 'cells')
    N = size(cellTemp)

    allocate(self % cells(N))
    self % cells % idx = cellTemp

    ! Load pointers and convert cell ID to IDX
    do i = 1, N
      ! Convert cell ID to IDX
      self % cells(i) % idx = cells % getIdx(self % cells(i) % idx)
      self % cells(i) % ptr => cells % getPtr(self % cells(i) % idx)
    end do

    ! Create fill array
    allocate(fill(N+2))

    do i = 1, N
      fill(i) = cells % getFill(self % cells(i) % idx)
    end do
    fill(N+1) = UNDEF_MAT
    fill(N+2) = OVERLAP_MAT

    ! Check for overlaps?
    call dict % getOrDefault(self % checkOverlap, 'checkOverlap', .false.)

    ! Create search indices to accelerate cell search
    allocate(self % searchIdx(N))
    self % searchIdx = [(i, i = 1, size(self % cells))]

    ! Initialise flag for blocking visit-ordered searches
    self % reordering = .false.

  end subroutine init

  !!
  !! Change the order in which cells are searched.
  !! Done by insertion sort.
  !! Should make cell searches more efficient, reordering to
  !! place the most frequently visited cells first.
  !! Insertion sort should be fast for nearly-sorted arrays.
  !!
  subroutine reorderCells(self)
    class(cellUniverse), intent(inout) :: self
    integer(shortInt)                  :: i, j, N, keyIdx
    integer(longInt)                   :: keyVisit

    associate(cells => self % cells, idx => self % searchIdx)
      N = size(cells)
      do i = 2, n
        keyIdx = idx(i)
        keyVisit = cells(keyIdx) % visits
        j = i - 1

        ! Shift elements right until correct position found
        do while (j >= 1 .and. cells(idx(j)) % visits < keyVisit)
          idx(j + 1) = idx(j)
          j = j - 1
        end do

        ! Insert element in correct position
        idx(j + 1) = keyIdx
      end do
    end associate

  end subroutine reorderCells

  !!
  !! Find local cell ID given a point
  !!
  !! See universe_inter for details.
  !!
  subroutine findCell(self, localID, cellIdx, r, u)
    class(cellUniverse), intent(inout)      :: self
    integer(shortInt), intent(out)          :: localID
    integer(shortInt), intent(out)          :: cellIdx
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    integer(shortInt)                       :: i
    integer(shortInt), dimension(:), allocatable :: array

    ! If being careful, only do exhaustive searches to check for overlap
    if (self % checkOverlap) then
      call self % findCellOverlap(localID, cellIdx, r, u)
      return
    end if

    ! Decide whether to use ordered array of cells or sorted array
    if (self % reordering) then
      array = [(i, i = 1, size(self % cells))]
    else
      array = self % searchIdx
    end if

    ! Search all cells
    do i = 1, size(self % cells)
      localID = array(i)
      if (self % cells(localID) % ptr % inside(r, u)) then
        cellIdx = self % cells(localID) % idx
        !$omp atomic
        self % cells(localID) % visits = self % cells(localID) % visits + 1
        
        if (mod(self % cells(localID) % visits, reorderFreq) == 0) then
          !$omp critical
          self % reordering = .true.
          call self % reorderCells()
          self % reordering = .false.
          !$omp end critical
        end if
        return

      end if
    end do

    ! If not found return undefined cell
    localID = i
    cellIdx = 0

  end subroutine findCell
  
  !!
  !! Find local cell ID given a point
  !! DEBUG procedure for detecting overlaps.
  !! Searches every cell and checks for more than one cell found.
  !! Significantly slower than alternatives in universes with
  !! many cells.
  !!
  !! See universe_inter for details.
  !!
  subroutine findCellOverlap(self, localID, cellIdx, r, u)
    class(cellUniverse), intent(inout)      :: self
    integer(shortInt), intent(out)          :: localID
    integer(shortInt), intent(out)          :: cellIdx
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    integer(shortInt)                       :: i, foundID, cellsFound
    character(100), parameter :: Here = 'findCellOverlap (cellUniverse_class.f90)'

    cellsFound = 0

    ! Search all cells
    do i = 1, size(self % cells)
      localID = i
      if (self % cells(localID) % ptr % inside(r, u)) then
        foundID = localID
        cellIdx = self % cells(localID) % idx
        cellsFound = cellsFound + 1
        !$omp atomic
        self % cells(localID) % visits = self % cells(localID) % visits + 1
        
      end if
    end do

    ! If not found return undefined cell
    if (cellsFound == 0) then
      localID = i
      cellIdx = 0
    ! If more than one found, return overlap cell
    elseif (cellsFound > 1) then
      localID = i + 1
      cellIdx = 0
    else
      localID = foundID
    end if

  end subroutine findCellOverlap
  
  !!
  !! Find local cell ID given a cell and its neighbour list
  !!
  subroutine findCellNeighb(self, localID, cellIdx, r, u, initID, foundNeighb)
    class(cellUniverse), intent(inout)      :: self
    integer(shortInt), intent(out)          :: localID
    integer(shortInt), intent(out)          :: cellIdx
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    integer(shortInt), intent(in)           :: initID
    logical(defBool), intent(out)           :: foundNeighb
    integer(shortInt)                       :: i

    ! If not found return undefined cell
    cellIdx = 0
    foundNeighb = .false.

    if (initID < 1) return

    associate(neighb => self % cells(initID) % neighb)

      ! Search neighbours if starting cell is given
      if (neighb % getSize() > 0) then
        do i = 1, neighb % getSize()

          localID = neighb % get(i)

          if (self % cells(localID) % ptr % inside(r, u)) then
            cellIdx = self % cells(localID) % idx
            foundNeighb = .true.
            !$omp atomic
            self % cells(localID) % visits = self % cells(localID) % visits + 1
            if (mod(self % cells(localID) % visits, reorderFreq) == 0) then
              !$omp critical
              self % reordering = .true.
              call self % reorderCells()
              self % reordering = .false.
              !$omp end critical
            end if
            return

          end if

        end do
      end if

    end associate

  end subroutine findCellNeighb

  !!
  !! Return distance to the next boundary between local cells in the universe
  !!
  !! See universe_inter for details.
  !!
  !! Errors:
  !!   fatalError if in UNDEFINED cell or OVERLAP cell
  !!
  subroutine distance(self, d, surfIdx, coords)
    class(cellUniverse), intent(inout) :: self
    real(defReal), intent(out)         :: d
    integer(shortInt), intent(out)     :: surfIdx
    type(coord), intent(in)            :: coords
    integer(shortInt)                  :: localID
    character(100), parameter :: Here = 'distance (cellUniverse_class.f90)'

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
    class(cellUniverse), intent(inout) :: self
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
    class(cellUniverse), intent(in) :: self
    type(coord), intent(in)         :: coords
    real(defReal), dimension(3)     :: offset

    ! There is no cell offset
    offset = ZERO

  end function cellOffset

  !!
  !! Return to uninitialised state
  !!
  subroutine kill(self)
    class(cellUniverse), intent(inout) :: self
    integer(shortInt)                  :: i

    ! Superclass
    call kill_super(self)

    ! Local
    if(allocated(self % cells)) then
      do i = 1, size(self % cells)
        call self % cells(i) % neighb % kill()
      end do
      deallocate(self % cells)
    end if
    if(allocated(self % searchIdx)) deallocate(self % searchIdx)

  end subroutine kill


end module cellUniverse_class
