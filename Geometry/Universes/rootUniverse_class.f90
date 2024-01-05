module rootUniverse_class

  use numPrecision
  use universalVariables, only : OUTSIDE_MAT
  use genericProcedures,  only : fatalError, numToChar
  use dictionary_class,   only : dictionary
  use coord_class,        only : coord
  use charMap_class,      only : charMap
  use surface_inter,      only : surface
  use surfaceShelf_class, only : surfaceShelf
  use cylinder_class,     only : cylinder
  use cell_inter,         only : cell
  use cellShelf_class,    only : cellShelf
  use universe_inter,     only : universe, kill_super => kill, charToFill

  implicit none
  private

  ! Parameters
  integer(shortInt), parameter :: INSIDE_ID = 1, OUTSIDE_ID = 2

  !!
  !! A top level (root) universe of geometry
  !!
  !! Is composed of two regions. Inside and outside separated by a single surface.
  !! Inside is the -ve halfspace of the boundary surface
  !! +ve halfspace is OUTSIDE
  !! Filling can be universe given by ID (`u<67` syntax) or a material given by name (e.g. 'fuel')
  !!
  !! Local ID 1 is inside. 2 is outside.
  !!
  !! Sample Input Dictionary:
  !!   root { type rootUniverse;
  !!          id 7;
  !!          border 78;   // Boundary surface
  !!          fill u<17>;  // Inside filling
  !!        }
  !!
  !!
  !! Public Members:
  !!   surf    -> Pointer to the boundary surface
  !!   surfIdx -> Index of the boundary surface
  !!
  !! Interface:
  !!   Universe interface
  !!   border -> Return surfIdx of the boundary surface
  !!
  type, public, extends(universe) :: rootUniverse
    class(surface), pointer :: surf => null()
    integer(shortInt)       :: surfIdx = 0
  contains
    ! Superclass procedures
    procedure :: init
    procedure :: kill
    procedure :: findCell
    procedure :: distance
    procedure :: cross
    procedure :: cellOffset

    ! Subclass procedures
    procedure :: border

  end type rootUniverse


contains

  !!
  !! Initialise Universe
  !!
  !! See universe_inter for details.
  !!
  !! Errors:
  !!   fatalError for invalid input
  !!
  subroutine init(self, fill, dict, cells, surfs, mats)
    class(rootUniverse), intent(inout)                        :: self
    integer(shortInt), dimension(:), allocatable, intent(out) :: fill
    class(dictionary), intent(in)                             :: dict
    type(cellShelf), intent(inout)                            :: cells
    type(surfaceShelf), intent(inout)                         :: surfs
    type(charMap), intent(in)                                 :: mats
    integer(shortInt)                                         :: id
    character(nameLen)                                        :: name
    character(100), parameter :: Here = 'init (rootUniverse_class.f90)'

    ! Setup the base class
    ! With: id, origin rotations...
    call self % setupBase(dict)

    ! Make sure root does not contain neither origin nor rotation
    if (dict % isPresent('origin')) then
      call fatalError(Here, 'Origin is not allowed. Centre of the root universe is &
                            &always (0.0 0.0 0.0).')
    else if (dict % isPresent('rotation')) then
      call fatalError(Here, 'Rotation is not allowed. Root universe cannot be rotated.')

    end if

    ! Get boundary surface
    call dict % get(id, 'border')

    if ( id == 0) then
      call fatalError(Here, 'Invalid border surface ID: 0')
    else if (id < 0) then
      call fatalError(Here, 'Border must be given as +ve ID. Inside is always in &
                           &-ve halfspace. Was given: '//numToChar(id))
   end if

   self % surfIdx = surfs % getIdx(id)
   self % surf => surfs % getPtr(self % surfIdx)

   ! Create fill array
   allocate(fill(2))
   fill(OUTSIDE_ID) = OUTSIDE_MAT
   call dict % get(name, 'fill')
   fill(INSIDE_ID) = charToFill(name, mats, Here)

  end subroutine init

  !!
  !! Find local cell ID given a point
  !!
  !! See universe_inter for details.
  !!
  subroutine findCell(self, localID, cellIdx, r, u)
    class(rootUniverse), intent(inout)      :: self
    integer(shortInt), intent(out)          :: localID
    integer(shortInt), intent(out)          :: cellIdx
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u

    cellIdx = 0
    ! Check halfspace
    if (self % surf % halfspace(r, u)) then
      localID = OUTSIDE_ID

    else
      localID = INSIDE_ID

    end if


  end subroutine findCell

  !!
  !! Return distance to the next boundary between local cells in the universe
  !!
  !! See universe_inter for details.
  !!
  subroutine distance(self, d, surfIdx, coords)
    class(rootUniverse), intent(inout) :: self
    real(defReal), intent(out)         :: d
    integer(shortInt), intent(out)     :: surfIdx
    type(coord), intent(in)            :: coords
    character(100), parameter :: Here = 'distance (rootUniverse_class.f90)'

    surfIdx = self % surfIdx
    d = self % surf % distance(coords % r, coords % dir)

  end subroutine distance

  !!
  !! Cross between local cells
  !!
  !! See universe_inter for details.
  !!
  !! Errors:
  !!   fatalError if surface from distance is not MOVING_IN or MOVING_OUT
  !!
  subroutine cross(self, coords, surfIdx)
    class(rootUniverse), intent(inout) :: self
    type(coord), intent(inout)         :: coords
    integer(shortInt), intent(in)      :: surfIdx
    character(100), parameter :: Here = 'cross (rootUniverse_class.f90)'

    ! Cross by cell finding in case of significant undershoots
    call self % findCell(coords % localID, coords % cellIdx, coords % r, coords % dir)

  end subroutine cross

  !!
  !! Return offset for the current cell
  !!
  !! See universe_inter for details.
  !!
  function cellOffset(self, coords) result (offset)
    class(rootUniverse), intent(in) :: self
    type(coord), intent(in)         :: coords
    real(defReal), dimension(3)     :: offset

    ! There is no cell offset
    offset = ZERO

  end function cellOffset

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(rootUniverse), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Kill local
    self % surf    => null()
    self % surfIdx = 0

  end subroutine kill

  !!
  !! Return Index of the boundary surface
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   Integer IDX of the boundary (border) surface
  !!
  pure function border(self) result(idx)
    class(rootUniverse), intent(in) :: self
    integer(shortInt)               :: idx

    idx = self % surfIdx

  end function border


end module rootUniverse_class
