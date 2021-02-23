module uniFills_class

  use numPrecision
  use universalVariables, only : targetNotFound, OUTSIDE_MAT
  use genericProcedures,  only : fatalError, numToChar, linFind
  use intMap_class,       only : intMap

  implicit none
  private

  ! Parameters
  integer(shortInt), parameter :: UNSET_ROOT = -8, USED = 1, UNUSED = 0

  !!
  !! Local Type to group information about a single universe
  !!
  !! Public Members:
  !!   id   -> ID of the universe
  !!   fill -> Array that maps materials and universe fills to each local universe cell (universe
  !!     fill is given as -uniIdx; material fill as matIdx )
  !!
  !! Note: During `uniFills` build universe filling may be stored as -uniID instead of IDX!
  !!
  type, public :: fillInfo
    integer(shortInt)                            :: id = 0
    integer(shortInt), dimension(:), allocatable :: fill
  end type fillInfo

  !!
  !! A Temporary structure to store universe fillings & check correctness
  !!
  !! Is used to transfer information about universe nesting structure decoupled from
  !! the spatial subdivision. Will be used to constrct geomGraph.
  !!
  !! TODO: The recursive procedures that do checks may not be the best from the point of
  !!   view of efficiency. Investigate if it is a problem.
  !!
  !! Public Mambers:
  !!  root -> Index of the root universe
  !!  uni  -> Array of `fillInfo` with data for each universe. Universe indices are the
  !!    same on the uni array on on the universe Shelf.
  !!
  !! Interface:
  !!   init          -> Set size to # of universes
  !!   addUniverse   -> Add universe entry
  !!   finishBuild   -> Translate uniIDs to uniIdxs and check that data for all universes was given
  !!   setRoot       -> Set Root universe index
  !!   hasCycles     -> Returns true if there is recursion in the definition
  !!   maxNesting    -> Return depth of nesting
  !!   nestedOutside -> Return true if outside material is present below root universe
  !!   unusedUniverses -> Return array of uniIdx that are defined but not used in geometry
  !!   kill          -> Return to uninitialised state
  !!
  type, public :: uniFills
    integer(shortInt)                         :: root = UNSET_ROOT
    type(fillInfo), dimension(:), allocatable :: uni
  contains
    ! Build procedures
    procedure :: init
    procedure :: addUniverse
    procedure :: finishBuild
    procedure :: setRoot
    procedure :: kill

    ! Structure checking
    procedure :: hasCycles
    procedure :: maxNesting
    procedure :: nestedOutside
    procedure :: unusedUniverses
    procedure :: countInstances

    ! Private procedures
    procedure, private :: isRepeated
    procedure, private :: countDepth
    procedure, private :: outsideBelow
    procedure, private :: collectUsed
    procedure, private :: countInstancesBelow

  end type uniFills

contains

  !!
  !! Set number of universes and allocate space for data
  !!
  !! Args:
  !!   N [in] -> +ve number of universes
  !!
  !! Errors:
  !!   fatalError if number of universes is -ve or 0
  !!
  subroutine init(self, N)
    class(uniFills), intent(inout) :: self
    integer(shortInt), intent(in)  :: N
    character(100), parameter :: Here = 'init (uniFills_class.f90)'

    if (N <= 0) call fatalError(Here, 'Given not +ve number of universes: '//numToChar(N))
    allocate(self % uni(N))

  end subroutine init

  !!
  !! Add a universe entry
  !!
  !! Args:
  !!   idx [in]     -> Universe index
  !!   id [in]      -> Universe ID
  !!   fill [inout] -> Allocatable array with fills. Will be deallocated on exit! (move_alloc)
  !!
  !! Errors:
  !!   fatalError if idx is outside the range
  !!
  subroutine addUniverse(self, idx, id, fill)
    class(uniFills), intent(inout)                              :: self
    integer(shortInt), intent(in)                               :: idx
    integer(shortInt), intent(in)                               :: id
    integer(shortInt), dimension(:), allocatable, intent(inout) :: fill
    integer(shortInt) :: top
    character(100), parameter :: Here = 'addUniverse (uniFills_class.f90)'

    ! Get top
    if (allocated(self % uni)) then
      top = size(self % uni)
    else
      top = 0
    end if

    ! Check index
    if (idx <= 0 .or. idx > top) then
      call fatalError(Here, 'Universe index: '//numToChar(idx)//' is outside the &
                            &valid range: 1->'//numToChar(top))
    end if

    ! Check for repetition
    if (allocated(self % uni(idx) % fill)) then
      call fatalError(Here, 'Cannot load universe (id '//numToChar(id)// ') with idx: '//&
                             numToChar(idx)//' becouse universe (id '//&
                             numToChar(self % uni(idx) % id)//') is alrady at that index.')
    end if

    ! Load data
    self % uni(idx) % id = id
    call move_alloc(fill, self % uni(idx) % fill)

  end subroutine addUniverse

  !!
  !! Finish build phase
  !!
  !! Change all uniIDs to uniIDx
  !! Check that all universe were loaded
  !!
  !! Args:
  !!   map [in] -> Map of uniIDs to uniIDXs
  !!
  !! Errors:
  !!   fatalError if some uniIdxs are empty (not all data was loaded)
  !!
  subroutine finishBuild(self, map)
    class(uniFills), intent(inout) :: self
    type(intMap), intent(in)       :: map
    integer(shortInt)              :: i, j, uniIdx
    character(100), parameter :: Here = 'finishBuild (uniFills_class.f90)'

    ! Check for gaps
    do i = 1, size(self % uni)
      if (.not.allocated(self % uni(i) % fill)) then
        call fatalError(Here, 'No data was loaded for universe with idx: '//numToChar(i))
      end if

    end do

    ! Translate IDs to IDXs
    do j = 1, size(self % uni) ! Loop over universes
      do i = 1, size(self % uni(j) % fill) ! Loop over local cells
        if (self % uni(j) % fill(i) < 0) then ! It is universe
          uniIdx = map % get(abs(self % uni(j) % fill(i)))
          self % uni(j) % fill(i) = -uniIdx

        end if

      end do
    end do

  end subroutine finishBuild

  !!
  !! Set root universe by IDX
  !!
  !! Args:
  !!   idx [in] -> uniIdx of the root universe (+ve integer)
  !!
  !! Errors:
  !!   fatalError if idx does not correspond to a defined universe
  !!
  subroutine setRoot(self, idx)
    class(uniFills), intent(inout) :: self
    integer(shortInt), intent(in)  :: idx
    integer(shortInt) :: top
    character(100), parameter :: Here = 'setRoot (uniFills_class.f90)'

    ! Get top
    if (allocated(self % uni)) then
      top = size(self % uni)
    else
      top = 0
    end if

    ! Check index
    if (idx <= 0 .or. idx > top) then
      call fatalError(Here, 'Root universe index: '//numToChar(idx)//' is outside the &
                            &valid range: 1->'//numToChar(top))
    end if

    ! Set value
    self % root = idx

  end subroutine setRoot

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(uniFills), intent(inout) :: self

    self % root = UNSET_ROOT
    if(allocated(self % uni)) deallocate(self % uni)

  end subroutine kill

  !!
  !! Returns true if there are cycles in geometry graph
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   True is there are cycles (recursion) in the definition
  !!
  !! Errors:
  !!   fatalError if root has not been set
  !!
  function hasCycles(self)
    class(uniFills), intent(in)     :: self
    logical(defBool)                :: hasCycles
    integer(shortInt), dimension(0) :: empty
    character(100), parameter :: Here = 'hasCycles (uniFills_class.f90)'

    ! Check that root is set
    if (self % root == UNSET_ROOT) then
      call fatalError(Here, 'Root universe has not been set.')
    end if

    hasCycles = self % isRepeated(empty, self % root)

  end function hasCycles

  !!
  !! Return maximum nesting level in the geometry structure
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   Maximum level of nesting
  !!
  !! Errors:
  !!   fatalError if root has not been set
  !!
  function maxNesting(self) result(N)
    class(uniFills), intent(in) :: self
    integer(shortInt)           :: N
    character(100), parameter :: Here = 'maxNesting (uniFills_class.f90)'

    ! Check that root is set
    if (self % root == UNSET_ROOT) then
      call fatalError(Here, 'Root universe has not been set.')
    end if

    N = self % countDepth(self % root)

  end function maxNesting

  !!
  !! Return true if outside appears below the root
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   True if outside appears below the root
  !!
  !! Errors:
  !!   fatalError if root is not set
  !!
  function nestedOutside(self) result(hasIt)
    class(uniFills), intent(in) :: self
    logical(defBool)            :: hasIt
    integer(shortInt)           :: i, fill
    character(100), parameter   :: Here = 'nestedOutside (uniFills_class.f90)'

    ! Check that root is set
    if (self % root == UNSET_ROOT) then
      call fatalError(Here, 'Root universe has not been set.')
    end if

    ! Serach nested universes
    hasIt = .false.
    do i = 1, size(self % uni(self % root) % fill)
      fill = self % uni(self % root) % fill(i)
      if (fill < 0) hasIt = self % outsideBelow(abs(fill))

      ! Return if found outside
      if (hasIt) return
    end do

  end function nestedOutside

  !!
  !! Return array of unused universes (by ID)
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   Array that contains unused universes (given as uniID)
  !!
  !! Errors:
  !!   fatalError is root is not set
  !!
  function unusedUniverses(self) result(list)
    class(uniFills), intent(in)                   :: self
    integer(shortInt), dimension(:), allocatable  :: list
    type(intMap)                                  :: usedUni
    logical(defBool), dimension(size(self % uni)) :: unusedMask
    integer(shortInt)                             :: i
    character(100), parameter   :: Here = 'unusedUniverses (uniFills_class.f90)'

    ! Check that root is set
    if (self % root == UNSET_ROOT) then
      call fatalError(Here, 'Root universe has not been set.')
    end if

    ! Collect used universes
    call usedUni % add(self % root, USED)
    call self % collectUsed(usedUni, self % root)

    ! Set logical mask
    unusedMask = .true.
    i = usedUni % begin()
    do while (i /= usedUni % end())
      unUsedMask(usedUni % atKey(i)) = .false.

      i = usedUni % next(i)
    end do

    ! Return unused
    list = pack(self % uni % id, unusedMask)

  end function unusedUniverses

  !!
  !! Get count of instances of diffrent universes in the geometry
  !!
  !! Args:
  !!   map [out] -> Map of uniIdx to number of instances. If uniIdx is not present in the
  !!     geometry, it is present in the map with value 0!
  !!
  subroutine countInstances(self, map)
    class(uniFills), intent(in) :: self
    type(intMap), intent(out)   :: map
    integer(shortInt)           :: i
    character(100), parameter :: Here = 'countInstances (uniFills_class.f90)'

    ! Check that root is set
    if (self % root == UNSET_ROOT) then
      call fatalError(Here, 'Root universe has not been set.')
    end if

    ! Initialise map with all uniIdx and 0 count
    do i = 1, size(self % uni)
      call map % add(i, 0)
    end do

    ! Perform count via recursion
    call self % countInstancesBelow(map, self % root)

  end subroutine countInstances

  !!
  !! Return true if there are repetitions in geometry graph
  !!
  !! Path to the current node is list of uniIdxs from root to the current node (including root)
  !! Current node is a universe given by idx.
  !!
  !! Args:
  !!   path [in] -> Universes contained in a path to the node given by idx
  !!   idx [in]  -> Current node universe given by idx
  !!
  !! Result:
  !!   True if current node or any universe in path is contained in the current universe
  !!   or below it.
  !!
  !! Error:
  !!   fatalError if idx is invalid
  !!
  recursive function isRepeated(self, path, idx) result(isIt)
    class(uniFills), intent(in)                 :: self
    integer(shortInt), dimension(:), intent(in) :: path
    integer(shortInt), intent(in)               :: idx
    logical(defBool)                            :: isIt
    integer(shortInt)                           :: i, pos, fill
    character(100), parameter :: Here = 'isRepeated (uniFills_class.f90)'

    ! Check if index is valid
    if (idx <= 0 .or. idx > size(self % uni)) then
      call fatalError(Here, 'Universe index: '//numToChar(idx)//' is invalid must &
                            &be in 1-'//numToChar(size(self % uni)))
    end if

    ! Loop over content
    isIt = .false.
    do i = 1, size(self % uni(idx) % fill)
      fill = self % uni(idx) % fill(i)

      if ( fill < 0) then
        ! Check if fill is the in path including current node
        pos = linFind(path, abs(fill))
        isIt = pos /= targetNotFound .or. abs(fill) == idx

        ! If repetition was not found search lower levels
        if (.not.isIt) isIt = self % isRepeated([path, idx], abs(fill))

        ! Return if repetition is found
        if (isIt) return
      end if

    end do

  end function isRepeated

  !!
  !! Return number of nesting levels below the node
  !!
  !! Args:
  !!   idx [in] -> Index of the current universe node
  !!
  !! Result:
  !!   Nesting levels below given universe including the given universe
  !!   (thus at minimum 1)
  !!
  !! Errors:
  !!   fatalError if index is invalid
  !!
  recursive function countDepth(self, idx) result(N)
    class(uniFills), intent(in)     :: self
    integer(shortInt), intent(in)   :: idx
    integer(shortInt)               :: N
    integer(shortInt)               :: i, fill
    character(100), parameter :: Here = 'countDepth (uniFills_class.f90)'

    ! Check if index is valid
    if (idx <= 0 .or. idx > size(self % uni)) then
      call fatalError(Here, 'Universe index: '//numToChar(idx)//' is invalid must &
                            &be in 1-'//numToChar(size(self % uni)))
    end if

    N = 1
    ! Loop over filling
    do i = 1, size(self % uni(idx) % fill)
      fill = self % uni(idx) % fill(i)
      ! If universe compare with current maximum depth
      if (fill < 0) N = max(N, 1 + self % countDepth(abs(fill)))

    end do

  end function countDepth

  !!
  !! Return true if outside is contained in the current universe or below it
  !!
  !! Args:
  !!   idx [in] -> Index of the current universe
  !!
  !! Result:
  !!   True if universe under idx contains outside or outisde is present below it
  !!
  !! Errors:
  !!   fatalError if idx is invalid
  !!
  recursive function outsideBelow(self, idx) result(hasIt)
    class(uniFills), intent(in)   :: self
    integer(shortInt), intent(in) :: idx
    logical(defBool)              :: hasIt
    integer(shortInt)             :: i, fill
    character(100), parameter :: Here = 'outsideBelow (uniFills_class.f90)'

    ! Check if index is valid
    if (idx <= 0 .or. idx > size(self % uni)) then
      call fatalError(Here, 'Universe index: '//numToChar(idx)//' is invalid must &
                            &be in 1-'//numToChar(size(self % uni)))
    end if

    ! Loop over local cells
    do i = 1, size(self % uni(idx) % fill)
      fill = self % uni(idx) % fill(i)
      hasIt = fill == OUTSIDE_MAT

      ! If fill is nested universe serach there
      if(fill < 0) hasIt = self % outsideBelow(abs(fill))

      ! Return if outside was found
      if(hasIt) return

    end do

  end function outsideBelow

  !!
  !! Add nested universes to the map (used as set)
  !!
  !! Args:
  !!   set [inout] -> Set of all used universes. Is a intMap (for now at least)
  !!   idx [in]    -> Index of the current universe
  !!
  !! Error:
  !!   fatalError if idx is invalid
  !!
  recursive subroutine collectUsed(self, set, idx)
    class(uniFills), intent(in)   :: self
    type(intMap), intent(inout)   :: set
    integer(shortInt), intent(in) :: idx
    integer(shortInt)             :: i, fill
    character(100), parameter :: Here = 'collectUsed (uniFills_class.f90)'

    ! Check if index is valid
    if (idx <= 0 .or. idx > size(self % uni)) then
      call fatalError(Here, 'Universe index: '//numToChar(idx)//' is invalid must &
                            &be in 1-'//numToChar(size(self % uni)))
    end if

    ! Loop over local cells
    do i = 1, size(self % uni(idx) % fill)
      fill = self % uni(idx) % fill(i)

      ! If fill is nested universe add it to set and go down the graph
      if(fill < 0) then
        call set % add(abs(fill), USED)
        call self % collectUsed(set, abs(fill))
      end if

    end do

  end subroutine collectUsed

  !!
  !! Count instances of diffrent univeres in the current universe or below it
  !!
  !! Args:
  !!   map [inout] -> Map of uniIdx to number of instances. Instances in the current universe
  !!     and below it are added to existing values
  !!   idx [in]    -> Index of the current universe
  !!
  recursive subroutine countInstancesBelow(self, map, idx)
    class(uniFills), intent(in)   :: self
    type(intMap), intent(inout)   :: map
    integer(shortInt), intent(in) :: idx
    integer(shortInt)             :: count, i, fill
    character(100), parameter :: Here = 'countInstancesBelow (uniFills_class.f90)'

    ! Check if index is valid
    if (idx <= 0 .or. idx > size(self % uni)) then
      call fatalError(Here, 'Universe index: '//numToChar(idx)//' is invalid must &
                            &be in 1-'//numToChar(size(self % uni)))
    end if

    ! Add this instance
    count = map % getOrDefault(idx, 0) + 1
    call map % add(idx, count)

    ! Loop over local cells and add nexted universes
    do i = 1, size(self % uni(idx) % fill)
      fill = self % uni(idx) % fill(i)

      ! If fill is nested universe count it and its contents
      if(fill < 0) then
        call self % countInstancesBelow(map, abs(fill))
      end if

    end do

  end subroutine countInstancesBelow

end module uniFills_class
