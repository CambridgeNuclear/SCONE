module csg_class

  use numPrecision
  use universalVariables
  use genericProcedures,    only : fatalError, numToChar, targetNotFound, linFind, printStart
  use dictionary_class,     only : dictionary
  use maps_class,           only : intMap
  use coord_class,          only : coord

  ! Material Names interface
  use nuclearData_inter,    only : nuclearData

  ! Surfaces
  use surfaceFactory_func,  only : new_surface
  use surface_inter,        only : surface, surfaceShelf

  ! Cells
  use cell_class,           only : cell, cellShelf, new_cell

  ! Universes
  use universeFactory_func, only : new_universe
  use universe_inter,       only : universe, universeShelf


  implicit none
  private


  !!
  !! Local helper wrapper around universe fill vector
  !! Vectors can have diffrent length so this allows us to store array of vectors
  !! +ve entries in fill vectors are cell IDs
  !! -ve entries in fill vector are uni IDs
  !!
  type, private :: uniFill
    integer(shortInt)                          :: uniID
    integer(shortInt)                          :: uniIdx
    integer(shortInt),dimension(:),allocatable :: fillV
  end type

  !!
  !! Simplified representation of CSG to perform checks and construct fillArray
  !!
  type, private :: csgTree
    type(uniFill), dimension(:), allocatable :: uniFills
    type(intMap)                             :: cellFillMap
    type(intMap)                             :: uniMap
  contains
    !! Procedures
    procedure :: init             => init_tree
    procedure :: isRecursive      => isRecursive_tree
    procedure :: maxNesting       => maxNesting_tree
    procedure :: cellCount        => cellCount_tree
    procedure :: hasNestedOutside => hasNestedOutside_tree
    procedure :: unusedUniverses  => unusedUniverses_tree

    !! Private recursion steps for procedures
    procedure :: appearsBelow     => appearsBelow_tree
    procedure :: nestingDepth     => nestingDepth_tree
    procedure :: countBelow       => countBelow_tree
    procedure :: containsOutside  => containsOutside_tree
    procedure :: setUsedUniverses => setUsedUniverses_tree
  end type csgTree

  !!
  !! Structure that allows to obtain fill for uniqueID and uniRootID for a transfer
  !! *** It is sub-optimal at the moment. There are gaps in id numbering
  !! *** Should be relativly easy to correct it
  !!
  type, public :: fillArray
    private
    integer(shortInt), dimension(:),allocatable :: fills
  contains
    procedure :: init    => init_fillArray
    procedure :: getFill => getFill_fillArray

    ! Private procedures
    procedure,private :: layoutUniverse => layoutUniverse_fillArray
  end type fillArray


  !!
  !! Representation of Constructive Solid Geometry (CSG)
  !! First entry on uShelf is a root Universe
  !! Will also store idx of boundary surface
  !!
  type, public :: csg
    integer(shortInt)   :: boundaryIdx = 0
    type(fillArray)     :: fills
    type(surfaceShelf)  :: sShelf
    type(cellShelf)     :: cShelf
    type(universeShelf) :: uShelf
  contains
    procedure :: init

    !! Build Macros
    procedure, private :: buildSurfaces
    procedure, private :: buildCells
    procedure, private :: buildUniverses
    procedure, private :: checkCellFills
    procedure, private :: findBoundarySurface
    procedure, private :: setBC
  end type csg

contains
  !!
  !! Build CSG representation from dictionary
  !!
  subroutine init(self,dict, materials)
    class(csg), intent(inout)                  :: self
    class(dictionary), intent(inout)           :: dict
    class(nuclearData), intent(in)             :: materials
    type(dictionary)                           :: surfDict
    type(dictionary)                           :: cellDict
    type(dictionary)                           :: uniDict
    type(intMap)                               :: fillMap
    type(uniFill),dimension(:),allocatable     :: uniFills
    integer(shortInt),dimension(:),allocatable :: BC
    type(csgTree)                              :: tree
    type(coord)          :: testco
    character(100),parameter :: Here = 'init (csg_class.f90)'
    integer(shortInt) :: cells, trans, nesting

    call printStart()

    print *, repeat('<>',50)
    ! *** START GEOMETRY BUILD
    print *, "/\/\ READING CSG GEOMETRY REPRESENTATION /\/\"

    ! Separate Dictionaries
    call dict % get(surfDict,'surfaces')
    call dict % get(cellDict,'cells')
    call dict % get(uniDict,'universes')

    ! Build defined surfaces
    call self % buildSurfaces(surfDict)

    ! Build defined cells
    call self % buildCells(fillMap, cellDict, materials)

    ! Build defined universes, put root universe under idx = 1
    call self % buildUniverses(uniFills, uniDict, materials)

    ! Remove extra memory from shelfs
    call self % sShelf % trimSize()
    call self % cShelf % trimSize()
    call self % uShelf % trimSize()

    ! Find boundary surface & set BC
    call self % findBoundarySurface(uniFills, fillMap)
    call self % setBC(dict)


    ! *** START GEOMETRY CHECKS
    print *, "CHECKING GEOMETRY:"

    ! Check that all universes filling cells are defined
    call self % checkCellFills(fillMap)

    ! Initialise csgTree
    call tree % init(fillMap, uniFills)

    ! Check for recursion
    if (tree % isRecursive()) then
      call fatalError(Here,'Recursion was found in geometry definition. &
                          & Universe cannot contain itself below itself.')
    else
      print *, "  Recursion in definition - NOT PRESENT!"
    end if

    ! Check maximum nesting levels
    nesting = tree % maxNesting()
    if (nesting > hardcoded_max_nest) then
      call fatalError(Here,'Nesting level: '// numToChar(nesting) //'> &
                          & max nesting'//numToChar(hardcoded_max_nest))
    else
      print *, "  Nesting level - FINE! "
    end if

    ! Check outside below root
    if (tree % hasNestedOutside()) then
      call fatalError(Here,'Cell with outside fill is present below root universe')
    else
      print *, "  Outside below root - NOT PRESENT!"
    end if

    ! *** START GEOMETRY INFORMATION
    print *, "GEOMETRY INFORMATION:"

    ! Read number of cells and transitions
    call tree % cellCount(cells,trans)

    ! Extract BC string
    call dict % get(BC,'boundary')

    ! Print geometry information
    print *,"  Nesting Levels: ",         numToChar(nesting)
    print *,"  Unique Cells: ",           numToChar(cells)
    print *,"  Unique Material Cells: ",  numToChar(cells - trans)
    print *,"  Nested Universes: ",       numToChar(trans)
    print *,"  Unused Universes by ID: ", numToChar(tree % unusedUniverses())
    print *,"  Boundary Surface ID: ",    numToChar(self % sShelf% shelf(self % boundaryIdx) % id())
    print *,"  Boundary Surface Type: ",   self % sShelf % shelf(self % boundaryIdx) % type()
    print *,"  Boundary Conditions: ",    numToChar(BC)
    ! TODO: ADD PROCEDURE FOR SURFACES TO PRINT BC AND PROVIDE BC INFORMATION

    ! *** BUILD FILL ARRAY
    print *, "BUILD FILL ARRAY:"
    call self % fills % init(tree)
    print *, "  Fill Array - DONE!"

    print *, "\/\/ FINISHED READING GEOMETRY \/\/"
    print *, repeat('<>',50)

  end subroutine init

  !!
  !! Build all surfaces in dictionary
  !!
  subroutine buildSurfaces(self,surfDict)
    class(csg), intent(inout)                   :: self
    class(dictionary), intent(inout)            :: surfDict
    character(nameLen),dimension(:),allocatable :: keys
    class(surface), allocatable                 :: tempSurf
    integer(shortInt)                           :: N, i, sIdx
    type(dictionary)                            :: tempDict

    ! Get all surfaces keys and number of surfaces
    call surfDict % keysDict(keys)
    N = size(keys)

    ! Initialise shelf
    call self % sShelf % init(N)

    ! Load surfaces on shelf
    print *, 'Building ', N, 'surfaces'

    do i=1,N
      ! Get dictionary and build surface from factory
      call surfDict % get(tempDict,keys(i))

      ! *** Quick fix allocate to temp allocatable
      allocate( tempSurf, source = new_surface(tempDict))
      call self % sShelf % addUnique( tempSurf, sIdx )

    end do

    print *, "DONE!"

  end subroutine buildSurfaces

  !!
  !! Build all cells in dictionary
  !!
  subroutine buildCells(self, fillMap ,cellDict, materials)
    class(csg), intent(inout)                   :: self
    type(intMap), intent(out)                   :: fillMap
    class(dictionary), intent(inout)            :: cellDict
    class(nuclearData), intent(in)              :: materials
    character(nameLen),dimension(:),allocatable :: keys
    type(cell)                                  :: tempCell
    type(dictionary)                            :: tempDict
    character(nameLen)                          :: filltype, matName
    integer(shortInt)                           :: N, i, id, fill, cellIdx, uniID
    character(100), parameter :: Here = 'buildCells (csg_class.f90)'

    ! Get all cell keys and number of cells
    call cellDict % keysDict(keys)
    N = size(keys)

    ! Initialise shelf and map
    call self % cShelf % init(N)
    call fillMap % init(N)
    print *, 'Building ', N, 'cells '

    ! Loop over all keys
    do i=1,N
      ! Get dictionary of a cell
      call cellDict % get(tempDict,keys(i))

      ! Build the cell
      tempCell = new_cell(tempDict, self % sShelf)
      call self % cShelf % add( tempCell, cellIdx )

      ! Get and verify cell id
      call tempDict % get(id,'id')
      if (id <= 0) call fatalError(Here,'Invalid cell id: '// numToChar(id))

      ! Determine fill of the cell and store appropertiate code on fillMap
      call tempDict % get(filltype,'filltype')

      select case(filltype)
        case('mat')
          ! Load material name
          call tempDict % get(matName,'mat')
          fill = materials % getIdx(matName)
          call fillMap % add(id, fill)

        case('outside')
          ! Load outside material
          fill = OUTSIDE_FILL
          call fillMap % add(id, fill)

        case('uni')
          ! Load nested universe ID
          call tempDict % get(uniID,'universe')
          if ( uniID < 1) call fatalError (Here,'Invalid fill universe ID' // numToChar(uniId))
          fill = -uniID
          call fillMap % add(id, fill)

        case default
          call fatalError(Here, 'Unrecognised fill type: ' // filltype )

      end select
    end do

    print *, "DONE!"

  end subroutine buildCells

  !!
  !! Build all universes in the dictionary
  !!
  subroutine buildUniverses(self, uniFills, uniDict, materials)
    class(csg), intent(inout)                          :: self
    type(uniFill),dimension(:),allocatable,intent(out) :: uniFills
    class(dictionary), intent(inout)                   :: uniDict
    class(nuclearData), intent(in)                     :: materials
    character(nameLen),dimension(:),allocatable        :: keys
    type(dictionary)                                   :: tempDict
    class(universe),allocatable                        :: tempUni
    integer(shortInt)                                  :: N, i, uniIdx, uniID

    ! Get all cell keys and number of cells
    call uniDict % keysDict(keys)
    N = size(keys)

    ! Initialise shelf and uniFills vector
    call self % uShelf % init(N)
    allocate(uniFills(N))
    print *, "Building ", N, "universes"

    ! Build root universe -I'm sorry. This is very ugly - MAK
    call uniDict % get(tempDict,'root')
    allocate( tempUni, source = new_universe(uniFills(1) % fillV,&
                                             tempDict,           &
                                             self % cShelf,      &
                                             self % sShelf,      &
                                             materials) )
    uniFills(1) % uniID = tempUni % id()
    call self % uShelf % add( tempUni, uniFills(1) % uniIdx)

    ! Remove root universe from the keys
    keys = pack(keys, keys /='root')

    ! Add the rest of universes
    do i=1,size(keys)
      call uniDict % get(tempDict, keys(i))
      allocate( tempUni, source = new_universe(uniFills(i+1) % fillV,&
                                               tempDict,           &
                                               self % cShelf,      &
                                               self % sShelf,      &
                                               materials) )
      uniFills(i+1) % uniID = tempUni % id()
      call self % uShelf % add( tempUni, uniFills(i+1) % uniIdx)

    end do

    print *, "DONE!"
  end subroutine buildUniverses

  !!
  !! Subroutine gives error if any cell is filled with an undefined universe
  !!
  subroutine checkCellFills(self,fillMap)
    class(csg), intent(in)   :: self
    type(intMap), intent(in) :: fillMap
    integer(shortInt)        :: i, cellID, fill, uniIdx
    character(100), parameter :: Here = 'checkCellFills (csg_class.f90)'

    do i=1,size(self % cShelf % shelf)
      ! Find cell Id
      cellID = self % cShelf % shelf(i) % id()
      ! Get fill
      fill = fillMap % get(cellID)
      if (fill < 0 ) then
        uniIdx = self % uShelf % getIdx(-fill)
        if (uniIdx == targetNotFound) then
          call fatalError(Here, 'Undefined uni: '//numToChar(-fill)//' in cell '//numToChar(cellID))

        end if
      end if
    end do
  end subroutine checkCellFills

  !!
  !! Set boundaryIdx to surface defining outside cell. Return error if outside is defined
  !! by multiple surfaces
  !!
  subroutine findBoundarySurface(self,uniFills, fillMap)
    class(csg), intent(inout)                  :: self
    type(uniFill),dimension(:),intent(in)      :: uniFills
    type(intMap),intent(in)                    :: fillMap
    integer(shortInt)                          :: N, i, N_out
    integer(shortInt)                          :: outsideCell, cellID, fill, cellIdx
    integer(shortInt),dimension(:),allocatable :: surfs
    character(100), parameter :: Here ='findBoundarySurface (csg_class.f90)'

    ! Find number of cells in root universe
    N = size( uniFills(1) % fillV)

    ! Set number of outside cells and outsideCell ID to 0
    N_out = 0
    outsideCell = 0

    ! Loop over cells in root universe
    do i=1,N
      fill = uniFills(1) % fillV(i)
      cellID = abs(fill)

      ! If fill is a cell, change fill to matIdx or universe index
      if (fill > 0) fill = fillMap % get(fill)

      ! If fill is outside set cellID of outside cell and increment outside cell counter
      if( fill == OUTSIDE_FILL ) then
        outsideCell = cellID
        N_out = N_out + 1

      end if
    end do

    ! Verify that single outside cell is present
    if (N_out < 1) then
      call fatalError(Here,'Outside cell is not present in root Universe')

    else if ( N_out > 1) then
      call fatalError(Here,'Only single outside cell is allowed for now')
    end if

    ! Find outside cell index and surfaces that define it
    cellIdx = self % cShelf % getIdx(outsideCell)
    surfs = self % cShelf % shelf(cellIdx) % getDef()

    ! Verify that only single surface is used for definition
    if( size(surfs) /= 1 ) call fatalError(Here,'Only single surface can define outside for now')

    ! Save index of boundary surface
    self % boundaryIdx = abs(surfs(1))

  end subroutine findBoundarySurface

  !!
  !! Read BC string from dictionary and apply it to boundary surface
  !! Boundary surface needs to be found before using this subroutine
  !!
  subroutine setBC(self, dict)
    class(csg), intent(inout)                    :: self
    class(dictionary), intent(in)                :: dict
    integer(shortInt), dimension(:), allocatable :: BC
    character(100),parameter :: Here ='setBC (csg_class.f90)'

    ! Verify that boundary surface index was set (contains valid idx)
    if(self % boundaryIdx < 1) then
      call fatalError(Here,'Boundary surface was not set before trying to load BC')
    end if

    ! Read boundary string
    call dict % get(BC,'boundary')

    ! Set boundary conditions on boundary surface
    associate ( bcSurf => self % sShelf % shelf( self % boundaryIdx))
      if( bcSurf % cannotBeBoundary()) then
        call fatalError(Here,'Boundary surface with ID: ' // numToChar(bcSurf % id()) // &
                             ' does not support Boundary Conditions')
      else
        call bcSurf % setBoundaryConditions(BC)

      end if
    end associate
  end subroutine setBC

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! csgTree functions
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Initialise csgTree from cellFill map, and uniFillVector array (uniFills)
  !! UniFills becomes unallocated upon exit
  !!
  subroutine init_tree(self, cellFillMap, uniFills)
    class(csgTree), intent(inout)                        :: self
    type(intMap), intent(inout)                          :: cellFillMap
    type(uniFill),dimension(:),allocatable,intent(inout) :: uniFills
    integer(shortInt)                                    :: i
    character(100), parameter :: Here = 'init_tree (csg_class.f90)'

    ! Initialise universe ID to Idx map
    call self % uniMap % init (size(uniFills))

    ! Check that uniIdx saved during universe build match with index
    do i=1,size(uniFills)
      if( uniFills(i) % uniIdx /= i) then
        call fatalError(Here,'uniIdx does not match with index on uniFills for uni ID: '&
                              // numToChar(uniFills(i) % uniID) )
      end if
      ! Load entry on universe ID to Idx map
      call self % uniMap % add(uniFills(i) % uniID, uniFills(i) % uniIdx)

    end do

    ! Load uniFills & cellFillMap
    self % cellFillMap = cellFillMap
    call move_alloc(uniFills, self % uniFills)

  end subroutine init_tree

  !!
  !! Returns true if there is a recursion in csg tree
  !! Recursion happens when universe is is contained within itself or beheth itself
  !!
  function isRecursive_tree(self) result(isIt)
    class(csgTree), intent(in)     :: self
    logical(defBool)               :: isIt
    integer(shortInt)              :: rootID
    integer(shortInt),dimension(0) :: emptyArray

    ! Find ID of root universe
    rootID = self % uniFills(1) % uniID

    ! Start recursive search
    isIt = self % appearsBelow(emptyArray, rootID)

  end function isRecursive_tree

  !!
  !! Function returns .true. if a any universe from uniIdsAbove is present in uniID or below
  !!
  recursive function appearsBelow_tree(self, uniIdsAbove, uniID) result(doesIt)
    class(csgTree), intent(in)                  :: self
    integer(shortInt), dimension(:), intent(in) :: uniIdsAbove
    integer(shortInt), intent(in)               :: uniID
    logical(defBool)                            :: doesIt
    integer(shortInt)                           :: N, i, uniIdx, fill, idx

    ! Find uniIdx of uniID
    uniIdx = self % uniMap % get(uniID)

    ! Find number of cells in uniID
    N = size(self % uniFills(uniIdx) % fillV)

    ! Loop over all cells - ignore outside fill
    do i=1,N
      fill = self % uniFills(uniIdx) % fillV(i)

      ! If fill is a cell change fill to matIdx or uniId
      if (fill > 0) fill = self % cellFillMap % get(fill)

      ! If fill is a universe look if it is present in uniIdsAbove
      ! It it is, set doesIt to .true.
      if (fill < 0) then
        idx = linFind(uniIdsAbove,abs(fill))
        doesIt = (idx /= targetNotFound)

        ! If repetition was not found search lower levels
        if(.not.doesIt) doesIt = self % appearsBelow([uniIdsAbove,uniId],abs(fill))

        ! Return if repetition was found
        if(doesIt) return

      end if
    end do

  end function appearsBelow_tree

  !!
  !! Calculates maximum nesting level of the CSG tree
  !!
  function maxNesting_tree(self) result(nest)
    class(csgTree), intent(in) :: self
    integer(shortInt)          :: nest
    integer(shortInt)          :: rootID

    ! Find ID of root universe
    rootID = self % uniFills(1) % uniID

    nest = self % nestingDepth(rootId)

  end function maxNesting_tree

  !!
  !! Recursion step for mexNesting_tree
  !! Adds one to maximum nestingDepth below the universe under uniID
  !!
  recursive function nestingDepth_tree(self,uniID) result (nest)
    class(csgTree), intent(in)    :: self
    integer(shortInt), intent(in) :: uniID
    integer(shortInt)             :: nest
    integer(shortInt)             :: N, i, uniIdx, fill, depth

    ! Find uniIdx of uniID
    uniIdx = self % uniMap % get(uniID)

    ! Find number of cells in uniID
    N = size(self % uniFills(uniIdx) % fillV)

    ! Loop over all cells - ignore outside fill
    depth = 0
    do i=1,N
      fill = self % uniFills(uniIdx) % fillV(i)

      ! If fill is a cell change fill to matIdx or uniId
      if (fill > 0) fill = self % cellFillMap % get(fill)

      ! If fill is a universe calculate its nesting depth and keep maximum
      if (fill < 0) then
        depth = max(depth, self % nestingDepth(abs(fill)))

      end if
    end do

    ! Add this universe to depth
    nest = depth + 1

  end function nestingDepth_tree

  !!
  !! Count number of cells and number of transfers from universe to universe in geometry
  !!
  subroutine cellCount_tree(self, cells, transfers)
    class(csgTree), intent(in)     :: self
    integer(shortInt), intent(out) :: cells
    integer(shortInt), intent(out) :: transfers
    integer(shortInt)              :: rootID

    ! Find ID of root universe
    rootID = self % uniFills(1) % uniID

    ! Set initial values to cells and transfersd
    cells     = 0
    transfers = 0

    ! Count cells & transfers in root and below
    call self % countBelow(cells, transfers, rootID)

  end subroutine cellCount_tree

  !!
  !! Adds cells and transfers in universe under uniID and below it to input cells & transfers
  !!
  recursive subroutine countBelow_tree(self, cells, transfers, uniID)
    class(csgTree), intent(in)       :: self
    integer(shortInt), intent(inout) :: cells
    integer(shortInt), intent(inout) :: transfers
    integer(shortInt), intent(in)    :: uniID
    integer(shortInt)                :: N, i, uniIdx, fill

    ! Find uniIdx of uniID
    uniIdx = self % uniMap % get(uniID)

    ! Find number of cells in uniID
    N = size(self % uniFills(uniIdx) % fillV)

    ! Increment number of cells
    cells = cells + N

    ! Loop over all cells - ignore outside fill
    do i=1,N
      fill = self % uniFills(uniIdx) % fillV(i)

      ! If fill is a cell change fill to matIdx or uniId
      if (fill > 0) fill = self % cellFillMap % get(fill)

      ! If fill is a universe add cells and tranfers below
      if (fill < 0) then
        ! Increment transfers count
        transfers = transfers + 1

        ! Count cells and transfers below
        call self % countBelow(cells, transfers, abs(fill))

      end if
    end do

  end subroutine countBelow_tree

  !!
  !! Returns true if there is an outside cell below root universe
  !!
  function hasNestedOutside_tree(self) result(doesIt)
    class(csgTree), intent(in) :: self
    logical(defBool)           :: doesIt
    integer(shortInt)          :: N, i, fill

    ! Find number of cells in root universe
    N = size(self % uniFills(1) % fillV)

    ! Initialise output to false
    doesIt = .false.

    ! Loop over all cells - ignore outside cells
    do i=1,N
      fill = self % uniFills(1) % fillV(i)

      ! If fill is a cell change fill to matIdx or uniId
      if (fill > 0) fill = self % cellFillMap % get(fill)

      ! If fill is a universe and contains outside below set doesIt = .true. and return
      if (fill < 0) then
        ! If there is outside below change fill to outside
        if( self % containsOutside(abs(fill))) then
          doesIt = .true.
          return

        end if
      end if
    end do

  end function hasNestedOutside_tree

  !!
  !! Returns true if universe under uniID contains outside in itself or below itself
  !!
  recursive function containsOutside_tree(self,uniId) result(doesIt)
    class(csgTree), intent(in)    :: self
    integer(shortInt), intent(in) :: uniId
    logical(defBool)              :: doesIt
    integer(shortInt)             :: N, i, uniIdx, fill

    ! Find uniIdx of uniID
    uniIdx = self % uniMap % get(uniID)

    ! Find number of cells in uniID
    N = size(self % uniFills(uniIdx) % fillV)

    ! Initialise output to false
    doesIt = .false.

    ! Loop over all cells
    do i=1,N
      fill = self % uniFills(uniIdx) % fillV(i)

      ! If fill is a cell change fill to matIdx or uniId
      if (fill > 0) fill = self % cellFillMap % get(fill)

      ! If fill is a universe and contains outside below change fill to outside
      if (fill < 0) then
        ! If there is outside below change fill to outside
        if( self % containsOutside(abs(fill))) then
          fill = 0
        end if
      end if

      ! If fill is outside set to true and return
      if (fill == 0 ) then ! is Outside fill
        doesIt = .true.
        return
      end if
    end do

  end function containsOutside_tree

  !!
  !! Function returns allocatable array that contains uniIds of unused universes
  !!
  function unusedUniverses_tree(self) result(unused)
    class(csgTree), intent(in)                        :: self
    integer(shortInt), dimension(:), allocatable      :: unused
    logical(defBool),dimension(size(self % uniFills)) :: uniVector
    integer(shortInt)                                 :: rootID

    ! Find ID of root universe
    rootID = self % uniFills(1) % uniID

    uniVector = .false.
    call self % setUsedUniverses(uniVector,rootId)

    ! Select unused universes
    unused = pack(self % uniFills % uniID, .not.uniVector)

  end function unusedUniverses_tree

  !!
  !! Given logical array coresponding in size to all avalible universes
  !! Set entry to .true. under index corresponding to universe present in uniId or below
  !!
  recursive subroutine setUsedUniverses_tree(self, uniVector, uniId)
    class(csgTree), intent(in)                    :: self
    logical(defBool), dimension(:), intent(inout) :: uniVector
    integer(shortInt), intent(in)                 :: uniId
    integer(shortInt)                             :: N, i, uniIdx, fill, fillIdx

    ! Find uniIdx of uniID
    uniIdx = self % uniMap % get(uniID)

    ! Set flag of uniIdx to true
    uniVector(uniIdx) = .true.

    ! Find number of cells in uniID
    N = size(self % uniFills(uniIdx) % fillV)

    ! Loop over all cells
    do i=1,N
      fill = self % uniFills(uniIdx) % fillV(i)

      ! If fill is a cell change fill to matIdx or uniId
      if (fill > 0) fill = self % cellFillMap % get(fill)

      ! If fill is a universe
      if (fill < 0) then
        ! Set corresponding entry in uniVEctor to .true.
        fillIdx = self % uniMap % get(abs(fill))
        uniVector(fillIdx) = .true.

        ! Check universes below
        call self % setUsedUniverses(uniVector,abs(fill))

      end if

    end do

  end subroutine setUsedUniverses_tree

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! fillArray procedures
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Build fillArray from csgTree
  !!
  subroutine init_fillArray(self,tree)
    class(fillArray), intent(inout)        :: self
    type(csgTree), intent(in)              :: tree
    integer(shortInt)                      :: cells, transfers
    type(uniFill),dimension(:),allocatable :: fillIdxs
    integer(shortInt)                      :: i,j, rootIdx, freeLoc
    character(100), parameter :: Here ='init_fillArray (csg_class.f90)'

    ! Recalculate number of cells and transfers to know required storage size
    ! Inefficient could be streamlined in a future
    call tree % cellCount(cells,transfers)

    ! Allocate storage
    allocate(self % fills(cells + transfers + 1))

    ! Copy fill vectors
    fillIdxs = tree % uniFills

    ! Loop over universes -> translate cell & uni Ids in fillIdxs into mat & uni IDXs
    do i=1,size(fillIdxs)
      ! Loop over all cells in a uiverse
      do j=1,size(fillIdxs(i) % fillV)
        associate (fill => fillIdxs(i) % fillV(j))

          ! Translate cell id to cell fill
          if( fill > 0 ) then ! fill is a cell
            fill = tree % cellFillMap % get(fill)
          end if

          ! If fill is universe change -uniID to -uniIdx
          if( fill < 0) then
            fill = -tree % uniMap % get(abs(fill))
          end if

        end associate
      end do
    end do

    ! Layout all universes
    rootIdx = 1
    freeLoc = 1
    call self % layoutUniverse(freeLoc,rootIdx,fillIdxs)

    ! Make a check for possible overflow when laying universe fills
    if ((freeLoc - 1) /= size(self % fills)) then
      call fatalError(Here,'Overflow of fillArray during filling')

    end if

  end subroutine init_fillArray

  !!
  !! Return a fill of a cell given by localID and uniRootID in coords
  !! If the cell is material fill is +ve and uniRootID is copied from coords
  !! If the cell is universe fill is -ve and uniRootID corresponds to the universe
  !!
  subroutine getFill_fillArray(self, coords, fill, uniRootID)
    class(fillArray), intent(in)   :: self
    type(coord), intent(in)        :: coords
    integer(shortInt), intent(out) :: fill
    integer(shortInt), intent(out) :: uniRootID
    integer(shortInt)              :: localID

    ! Read fill from fill array
    localID   = coords % localID
    uniRootID = coords % uniRootID

    fill = self % fills(uniRootID + localID)

    ! If fill is universe read uniRootID and read fill uniIdx
    if (fill < 0 ) then
      uniRootID = - fill
      fill      = - self % fills(uniRootID)
    end if

  end subroutine getFill_fillArray


  !!
  !! Set entries starting with freeLoc to fill of universe under uniIdx
  !! Use provided fillIdxs with contains fill of universes in terms of mat & uni idxs
  !!
  recursive subroutine layoutUniverse_fillArray(self, freeLoc, uniIdx, fillIdxs)
    class(fillArray), intent(inout)        :: self
    integer(shortInt), intent(inout)       :: freeLoc
    integer(shortInt), intent(in)          :: uniIdx
    type(uniFill),dimension(:), intent(in) :: fillIdxs
    integer(shortInt)                      :: N, startLoc, i, fill

    ! Find number of entries in universe unde uniIDx
    N = size (fillIdxs(uniIdx) % fillV)

    ! Save initial freeLocation
    startLoc = freeLoc

    ! Layout universe fill
    self % fills( startLoc) = uniIdx
    self % fills( startLoc + 1 : startLoc + N ) = fillIdxs(uniIdx) % fillV

    ! Update free Loc
    freeLoc = startLoc + N + 1

    ! Loop over layout values and layout nested universes
    do i = startLoc+1, startLoc + N
      fill = self % fills(i)

      if (fill < 0) then
        ! Save pointer to entry for the nested universe
        self % fills(i) = -freeLoc
        ! Layout nested universe
        call self % layoutUniverse(freeLoc, abs(fill), fillIdxs)

      end if
    end do

  end subroutine layoutUniverse_fillArray

end module csg_class
