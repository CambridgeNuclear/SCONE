module csg_class

  use numPrecision
  use universalVariables,  only : HARDCODED_MAX_NEST
  use display_func,        only : statusMsg, printSectionStart, printSectionEnd, printSeparatorLine
  use genericProcedures,   only : fatalError, numToChar
  use dictionary_class,    only : dictionary
  use charMap_class,       only : charMap
  use surfaceShelf_class,  only : surfaceShelf
  use surface_inter,       only : surface
  use cellShelf_class,     only : cellShelf
  use universeShelf_class, only : universeShelf
  use universe_inter,      only : universe
  use rootUniverse_class,  only : rootUniverse
  use uniFills_class,      only : uniFills
  use geomGraph_class,     only : geomGraph

  implicit none

  !!
  !! Representation of geometry spatial information
  !!
  !! Not exactly CSG in the standard sense, but name is short and convenient.
  !! It is used to store & initialise all geometry components.
  !!
  !! Sample Input Dictionary:
  !!   geometry {
  !!     #root 7#  // Optional universe ID of the universe
  !!     surfaces { < Surface Shelf definition>}
  !!     cells {<Cell Shelf definition>}
  !!     universes {<Universe SHelf definition>}
  !!     graph {<Geometry graphs definition>}
  !!     boundary (1 0 0 0 0 2);
  !!   }
  !!
  !! If `root` is not set a universe with name 'root' will be selected as a root universe.
  !!
  !! Public Members:
  !!   surfs     -> Shelf with all surfaces
  !!   cells     -> Shelf with all cells
  !!   unis      -> Shelf with all universes
  !!   graph     -> Geometry structure graph
  !!   borderIdx -> Index of the boundary surface
  !!   rootIdx   -> Index of the root universe
  !!
  !! Interface:
  !!   init -> Load geometry
  !!   kill -> Return to uninitialised state
  !!
  type, public :: csg
    type(surfaceShelf)  :: surfs
    type(cellShelf)     :: cells
    type(universeShelf) :: unis
    type(geomGraph)     :: graph
    integer(shortInt)   :: borderIdx = 0
    integer(shortInt)   :: rootIdx   = 0
  contains
    procedure :: init
    procedure :: kill
  end type csg


contains

  !!
  !! Initialise geometry representation
  !!
  !! Args:
  !!   sict [in] -> Dictionary with geometry definition
  !!   mats [in] -> Map of material names to matIdx
  !!   silent [in] -> Optional. Set to .true. to surpress console output (default .false.)
  !!
  !! Errors:
  !!   fatalError if geometry definition is invalid
  !!
  subroutine init(self, dict, mats, silent)
    class(csg), intent(inout)                    :: self
    class(dictionary), intent(in)                :: dict
    type(charMap), intent(in)                    :: mats
    logical(defBool), optional, intent(in)       :: silent
    logical(defBool)                             :: loud
    type(uniFills)                               :: fills
    integer(shortInt)                            :: rootId, nesting
    type(dictionary), pointer                    :: tempDict
    class(universe), pointer                     :: uni_ptr
    class(surface), pointer                      :: surf_ptr
    integer(shortInt), dimension(:), allocatable :: BC
    character(100), parameter :: Here = 'init (csg_class.f90)'

    ! Choose whether to display messages
    if (present(silent)) then
      loud = .not.silent
    else
      loud = .true.
    end if

    ! Print beggining
    if (loud) then
      call printSeparatorLine()
      call printSectionStart("READING GEOMETRY")
    end if

    ! Build Surfaces
    if (loud) call statusMsg("Building Surfaces")
    call self % surfs % init(dict % getDictPtr('surfaces'))
    if (loud) call statusMsg("DONE!")

    ! Build Cells
    if (loud) call statusMsg("Building Cells")
    call self % cells % init(dict % getDictPtr('cells'), self % surfs, mats)
    if (loud) call statusMsg("DONE!")

    ! Build Universes
    if(loud) call statusMsg("Building Universes")
    call self % unis % init(fills, &
                            dict % getDictPtr('universes'),&
                            self % cells, &
                            self % surfs, &
                            mats)
    if (loud) call statusMsg("DONE!")

    ! Select Root universe
    if (dict % isPresent('root')) then
      call dict % get(rootId, 'root')

    else
      tempDict => dict % getDictPtr('universes')
      tempDict => tempDict % getDictPtr('root')
      call tempDict % get(rootID, 'id')

    end if

    ! Set Root universe idx
    self % rootIdx = self % unis % getIdx(rootID)
    call fills % setRoot(self % rootIdx)

    ! Check that root universe is `rootUniverse`
    ! If it is set borderIdx
    uni_ptr => self % unis % getPtr(self % rootIdx)

    select type(uni_ptr)
    class is (rootUniverse)
      ! Get boundary surface
      self % borderIdx = uni_ptr % border()

    class default
      call fatalError(Here, 'Root universe with ID: '//numToChar(rootID)//' is not type &
                           &`rootUniverse`. Thus it cannot be the root of the geometry')
    end select

    ! Set boundary conditions
    call dict % get(BC, 'boundary')
    surf_ptr => self % surfs % getPtr(self % borderIdx)
    call surf_ptr % setBC(BC)

    ! Check validity of geometry structure
    if (loud) call statusMsg("CHECKING GEOMETRY:")
    ! Check for recursion
    if (fills % hasCycles()) then
      call fatalError(Here ,'There is recursion in the geometry nesting. &
                            &Universe cannot contain itself below itself.')
    else if (loud) then
      call statusMsg("  Recursion in definition - NOT PRESENT!")
    end if

    ! Check maximum nesting
    nesting = fills % maxNesting()
    if (nesting > HARDCODED_MAX_NEST) then
      call fatalError(Here,'Nesting level: '// numToChar(nesting) //'> &
                          & max nesting'//numToChar(HARDCODED_MAX_NEST))
    else if (loud) then
      call statusMsg("  Nesting level - FINE!")
    end if

    ! Check outside below root
    if (fills % nestedOutside()) then
      call fatalError(Here,'Cell with outside fill is present below root universe')
    else if (loud) then
      call statusMsg("  Outside below root - NOT PRESENT!")
    end if

    ! Build geometry Graph
    if (loud) call statusMsg("BUILDING GEOMETRY GRAPH")
    call self % graph % init(fills, dict % getDictPtr('graph'))
    if (loud) call statusMsg("DONE!")

    ! Print geometry information
    if (loud) then
      call statusMsg("GEOMETRY INFORMATION")
      call statusMsg("  Number of Surfaces: " // numToChar(self % surfs % getSize()))
      call statusMsg("  Number of Cells: " // numToChar(self % cells % getSize()))
      call statusMsg("  Number of Universes: " // numToChar(self % unis % getSize()))
      call statusMsg("  Nesting Levels: " // numToChar(nesting))
      call statusMsg("  Unique Cells: " // numToChar(self % graph % uniqueCells))
      call statusMsg("  Unused universes (ID): " // numToChar(fills % unusedUniverses()))
      call statusMsg("  Boundary Surface ID: " // numToChar(surf_ptr % id()))
      call statusMsg("  Boundary Surface Type: " // surf_ptr % myType())
      call statusMsg("  Boundary Conditions: " // numToChar(BC))
    end if

    ! Print End
    if (loud) then
      call printSectionEnd("FINISHED READING GEOMETRY")
      call printSeparatorLine
    end if

  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(csg), intent(inout) :: self

    ! Clean content
    call self % surfs % kill()
    call self % cells % kill()
    call self % unis % kill()
    call self % graph % kill()
    self % borderIdx = 0
    self % rootIdx = 0

  end subroutine kill

end module csg_class
