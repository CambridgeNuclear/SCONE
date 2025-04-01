module cellUniverse_class

  use numPrecision
  use universalVariables, only : UNDEF_MAT, NUDGE
  use genericProcedures,  only : fatalError, numToChar
  use dictionary_class,   only : dictionary
  use coord_class,        only : coord
  use charMap_class,      only : charMap
  use surfaceShelf_class, only : surfaceShelf
  use cell_inter,         only : cell
  use cellShelf_class,    only : cellShelf
  use universe_inter,     only : universe, kill_super => kill
  use omp_lib,            only : omp_test_lock, omp_unset_lock, &
                                 omp_init_lock, omp_destroy_lock, &
                                 omp_lock_kind
  implicit none
  private

  !!
  !! Local helper class to group cell data
  !!
  !! Public Members:
  !!   idx -> cellIdx of the cell in cellShelf
  !!   ptr -> Pointer to the cell
  !!
  type, private :: localCell
    integer(shortInt)    :: idx = 0
    integer(longInt)     :: visits = 0
    integer(shortInt), dimension(:), allocatable :: neighb
    class(cell), pointer :: ptr => null()
  end type localCell

  !!
  !! Representation of a universe via cells
  !!
  !! Each local cell in the universe corresponds to a cell given by an ID.
  !! An extra local cell is always defined inside the cellUniverse with UNDEF_MAT
  !! (undefined material) filling. If position is not in any user-defined cell, it is in this
  !! extra cell. Extra cell exists to enable plotting of geometry without fatalErrors.
  !!
  !! Sample Input Dictionary:
  !!   uni { type cellUniverse;
  !!         id 7;
  !!         # origin (2.0 0.0 0.0);    #
  !!         # rotation (23.0 0.0 0.0); #
  !!         cells ( 1 3 4);         }
  !!
  !! Note:
  !!   - Local IDs are assigned in order as in definition. In the example above local id would map
  !!   to following cell  [localID: cellID] {1: 1, 2: 3, 3: 4, 4: UNDEF }
  !!   - Cell overlaps are forbidden, but there is no check to find overlaps.
  !!
  !! Public Members:
  !!   cells -> Structure that stores cellIdx and pointers to the cells
  !!
  !! Interface:
  !!   universe interface
  !!
  type, public, extends(universe) :: cellUniverse
    type(localCell), dimension(:), allocatable :: cells
#ifdef _OPENMP
    integer(kind=omp_lock_kind), dimension(:), allocatable :: cellLocks
#endif
  contains
    ! Superclass procedures
    procedure :: init
    procedure :: kill
    procedure :: findCell
    procedure :: distance
    procedure :: cross
    procedure :: cellOffset
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
    allocate(fill(N+1))

    do i = 1, N
      fill(i) = cells % getFill(self % cells(i) % idx)
    end do
    fill(N+1) = UNDEF_MAT

    ! Create OMP locks for the cells
#ifdef _OPENMP
    allocate(self % cellLocks(N))
    do i = 1, N
      call OMP_init_lock(self % cellLocks)
    end do
#endif

  end subroutine init

  !!
  !! Find local cell ID given a point
  !!
  !! See universe_inter for details.
  !!
  subroutine findCell(self, localID, cellIdx, r, u, initID, newNeighb)
    class(cellUniverse), intent(inout)      :: self
    integer(shortInt), intent(out)          :: localID
    integer(shortInt), intent(out)          :: cellIdx
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    logical(defReal), intent(out), optional :: newNeighb
    integer(shortInt)                       :: initID, i

    newNeighb = .true.

    ! Search neighbours if starting cell is given
    if (present(initID)) then
            
      if ((initID > 0) .and. allocated(self % cells(initID) % neighb)) then
        do i = 1, size(self % cells(initID) % neighb)
        
          localID = self % cells(initID) % neighb(i)
        
          if (self % cells(localID) % ptr % inside(r, u)) then
            cellIdx = self % cells(localID) % idx
            newNeighb = .false.
            !$omp atomic
            self % cells(localID) % visits = self % cells(localID) % visits + 1
            return

          end if

        end do
      end if
    end if

    ! Search all cells
    do localID = 1, size(self % cells)
      if (self % cells(localID) % ptr % inside(r, u)) then
        cellIdx = self % cells(localID) % idx
        !$omp atomic
        self % cells(localID) % visits = self % cells(localID) % visits + 1
        return

      end if
    end do

    ! If not found return undefined cell
    ! Already set to localID (== size(self % cells) + 1) by the do loop
    cellIdx = 0

  end subroutine findCell

  !!
  !! Return distance to the next boundary between local cells in the universe
  !!
  !! See universe_inter for details.
  !!
  !! Errors:
  !!   fatalError if in UNDEFINED cell
  !!
  subroutine distance(self, d, surfIdx, coords)
    class(cellUniverse), intent(inout) :: self
    real(defReal), intent(out)         :: d
    integer(shortInt), intent(out)     :: surfIdx
    type(coord), intent(in)            :: coords
    integer(shortInt)                  :: localID
    character(100), parameter :: Here = 'distance (cellUniverse_class.f90)'

    localID = coords % localID

    if (localID > size(self % cells)) then
      call fatalError(Here, 'Particle is in undefined local cell. Local ID: '//numToChar(localID))
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
    logical(defBool)                   :: lockSet

    ! Keep initial cell ID to add to neighbour lists
    local0 = coords % localID

    ! NUDGE position slightly forward to escape surface tolerance
    ! and avoid calculating normal and extra dot-products
    coords % r = coords % r + coords % dir * NUDGE

    ! Find cell
    ! TODO: Some cell neighbout list
    call self % findCell(coords % localID, &
                         coords % cellIdx, &
                         coords % r,       &
                         coords % dir,     &
                         initID = local0,  &
                         newNeighb = newNeighb)

    ! Add to neighbour list
    if (newNeighb) then
      lockSet = .true.

#ifdef _OPENMP
      lockSet = OMP_test_lock(self % cellLocks(local0))
#endif
      if (lockSet) then
         call self % addNeighbour(local0, coords % localID)
#ifdef _OPENMP
         call OMP_unset_lock(self % cellLocks(local0))
#endif
      end if

#ifdef _OPENMP
      lockSet = OMP_test_lock(self % cellLocks(coords % localID))
#endif
      if (lockSet) then
         call self % addNeighbour(coords % localID, local0)
#ifdef _OPENMP
         call OMP_unset_lock(self % cellLocks(coords % localID))
#endif
      end if

    end if

  end subroutine cross

  !!
  !! Add one cell to another's neighbour list
  !!
  subroutine addNeighbour(self, id1, id2)
    class(cellUniverse), intent(inout)           :: self
    integer(shortInt), intent(in)                :: id1
    integer(shortInt), intent(in)                :: id2
    integer(shortInt)                            :: sz
    integer(shortInt), dimension(:), allocatable :: temp

    associate(neighb => self % cells(id1) % neighb)
    
      ! Size list
      sz = size(neighb)

      ! Handle case where array is not allocated
      if (sz1 == 0) then
        allocate(neighb)
        neighb(1) = id2
      else
      
        ! Ensure neighbour isn't already present
        if (any(neighb == id2)) return
        
        temp = neighb
        deallocate(neighb)
        allocate(neighb(sz + 1))
        neighb(sz + 1) = temp
      end if

    end associate

  end subroutine addNeighbour

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
  elemental subroutine kill(self)
    class(cellUniverse), intent(inout) :: self

    ! SUperclass
    call kill_super(self)

    ! Local
    if(allocated(self % cells)) deallocate(self % cells)
#ifdef _OPENMP
    if(allocated(self % cellLocks)) deallocate(self % cellLocks)
#endif

  end subroutine kill


end module cellUniverse_class
