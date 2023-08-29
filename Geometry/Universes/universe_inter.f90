module universe_inter

  use numPrecision
  use genericProcedures,  only : fatalError, numToChar, rotationMatrix, charToInt
  use dictionary_class,   only : dictionary
  use coord_class,        only : coord
  use charMap_class,      only : charMap
  use surfaceShelf_class, only : surfaceShelf
  use cellShelf_class,    only : cellShelf

  implicit none
  private


  ! Extandable methods
  public :: kill

  ! Universe utility functions
  public :: charToFill

  !!
  !! Abstract interface for all universes
  !!
  !! Universe represents a subdivision of the entire space into local cells.
  !!
  !! Universe can be associated with:
  !!   translation (to a new origin)
  !!   rotation (by Euler angles using ZXZ convention)
  !!
  !! Rotation is applied before translation to the origin.
  !!
  !! Sample Dictionary Input:
  !!   uni {
  !!     id 7;
  !!     #origin (1.0 0.0 0.1);#
  !!     #rotation (30.0 0.0 0.0);#
  !!     <subclass specific data>
  !!   }
  !!
  !! Private Members:
  !!   uniId   -> Id of the universe
  !!   uniIdx  -> Index of the universe
  !!   origin  -> Location of the origin of the univere co-ordinates in the frame of higher universe
  !!   rotMat  -> Rotation matrix for rotation with respect to the higher universe
  !!   rot     -> rotation flag. True is universe is rotated
  !!
  !! Interface:
  !!   id           -> Get Id of the universe
  !!   setId        -> Set Id of the universe
  !!   setIdx       -> Set index of the universe
  !!   setTransfrom -> Set origin and/or rotation of the universe. Rotation angles are in [deg]
  !!   init         -> Initialise universe and return fillArray with content of local cells.
  !!     Requires surface/cell shelfs and map of material names to matIdxs
  !!   kill         -> Return to uninitialised state
  !!   enter        -> Generate new coord after entering a universe from higher level coord.
  !!   findCell     -> Return local cell ID and cellIdx in cellShelf for a given position.
  !!   distance     -> Calculate the distance to the point where localID will change
  !!   cross        -> Assuming point is at the boundary between local cells, find next local cell.
  !!   cellOffset   -> Given coords with set localID, return cellOffset for that cell.
  !!
  type, public, abstract :: universe
    private
    integer(shortInt)             :: uniId  = 0
    integer(shortInt)             :: uniIdx = 0
    real(defReal), dimension(3)   :: origin = ZERO
    real(defReal), dimension(3,3) :: rotMat = ZERO
    logical(defBool)              :: rot    = .false.
  contains
    ! Build procedures
    procedure, non_overridable :: id
    procedure, non_overridable :: setId
    procedure, non_overridable :: setIdx
    procedure, non_overridable :: setTransform
    procedure(init), deferred  :: init
    procedure                  :: setupBase
    procedure                  :: kill

    ! Runtime procedures
    procedure, non_overridable      :: enter
    procedure(findCell), deferred   :: findCell
    procedure(distance), deferred   :: distance
    procedure(cross), deferred      :: cross
    procedure(cellOffset), deferred :: cellOffset
  end type universe

  abstract interface

    !!
    !! Initialise Universe
    !!
    !! Must ruturn a fill array, that contains content in each local cell of the universe.
    !! Array is indexed by local cell ID. Universe content is -uniID and material content
    !! is +ve matIdx.
    !!
    !! Args:
    !!   fill [out]    -> Fill array. Each position gives content for each local cell. Negative
    !!     values are universe IDs (uniID). Positive are material indexes (matIdx)
    !!   dict [in]     -> Dictionary with the universe definition
    !!   cells [inout] -> Shelf with all user-defined cells
    !!   surfs [inout] -> Shelf with all user-defined surfaces
    !!   mats [in]     -> Map of material names to corresponding matIdx
    !!
    subroutine init(self, fill, dict, cells, surfs, mats)
      import :: universe, shortInt, dictionary, &
                cellShelf, surfaceShelf, charMap
      class(universe), intent(inout)                            :: self
      integer(shortInt), dimension(:), allocatable, intent(out) :: fill
      class(dictionary), intent(in)                             :: dict
      type(cellShelf), intent(inout)                            :: cells
      type(surfaceShelf), intent(inout)                         :: surfs
      type(charMap), intent(in)                                 :: mats
    end subroutine init

    !!
    !! Find local cell ID given a point
    !!
    !! Given position and direction return localID for a cell and, if the cell is
    !! also defined on a cellShelf, return its index (cellIdx). If it isn't set
    !! cellIdx to 0.
    !!
    !! Args:
    !!   localID [out] -> Local ID for the given point
    !!   cellIdx [out] -> cellIdx in cellShelf, if the cell point is in is defined there.
    !!     If the cell exists only in the universe return 0.
    !!   r [in]        -> Position of a point
    !!   u [in]        -> Normalised direaction (norm2(u) = 1.0)
    !!
    !! Note: Self is intent(inout), but if a state of the universe is to be changed
    !!   it is necessary to consider issues related to parallel calculations with shared
    !!   memory.
    !!
    subroutine findCell(self, localID, cellIdx, r, u)
      import :: universe, shortInt, defReal
      class(universe), intent(inout)          :: self
      integer(shortInt), intent(out)          :: localID
      integer(shortInt), intent(out)          :: cellIdx
      real(defReal), dimension(3), intent(in) :: r
      real(defReal), dimension(3), intent(in) :: u
    end subroutine findCell

    !!
    !! Return distance to the next boundary between local cells in the universe
    !!
    !! In addition to distance surfIdx of the surface that will be crossed after a move by
    !! the distance is also returned. If surfIdx > 0 it means that the surface is defined
    !! on surfaceShelf with the index. If surfIdx < 0 it means that the surface is local to the
    !! universe.
    !!
    !! The returned surfIdx is a hint for `cross` procedure, that may allow to speed up transition
    !! processing.
    !!
    !! Args:
    !!   d [out]       -> Distance to the next surface
    !!   surfIdx [out] -> Index of the surface that will be crossed. If +ve than surface
    !!     is defined on surfaceSHelf. If -ve surface is local to this univerese.
    !!   coords [in]   -> Coordinates of the point inside the universe (after transformations
    !!     and with localID already set)
    !!
    !! Note: Self is intent(inout), but if a state of the universe is to be changed
    !!   it is necessary to consider issues related to parallel calculations with shared
    !!   memory.
    !!
    subroutine distance(self, d, surfIdx, coords)
      import :: universe, defReal, shortInt, coord
      class(universe), intent(inout)          :: self
      real(defReal), intent(out)              :: d
      integer(shortInt), intent(out)          :: surfIdx
      type(coord), intent(in)                 :: coords
    end subroutine distance

    !!
    !! Cross between local cells
    !!
    !! Procedure assumes that the point is ON THE SURFACE between cells within under/overshoot as
    !! a result of finate FP precision.
    !!
    !! Args:
    !!   coords [inout] -> Coordinates placed in the universe (after transformations and with
    !!     local ID set). On exit localID will be changed
    !!   surfIdx [in]   -> surfIdx from distance procedure, which hints which surface is beeing
    !!     crossed.
    !!
    !! Note: Self is intent(inout), but if a state of the universe is to be changed
    !!   it is necessary to consider issues related to parallel calculations with shared
    !!   memory.
    !!
    subroutine cross(self, coords, surfIdx)
      import :: universe, coord, shortInt
      class(universe), intent(inout) :: self
      type(coord), intent(inout)     :: coords
      integer(shortInt), intent(in)  :: surfIdx
    end subroutine cross

    !!
    !! Return offset for the current cell
    !!
    !! Args:
    !!   coords [in] -> Coordinates placed in the universe (after transformations and with
    !!     local ID set).
    !!
    !! Result:
    !!   Cell offset (3D position vector). Offset is applied before entering a nested universe
    !!   inside a local cell.
    !!
    function cellOffset(self, coords) result (offset)
      import :: universe, coord, defReal
      class(universe), intent(in) :: self
      type(coord), intent(in)     :: coords
      real(defReal), dimension(3) :: offset
    end function cellOffset

  end interface

contains

  !!
  !! Set universe ID
  !!
  !! Args:
  !!   id [in] -> Universe ID (id > 0)
  !!
  !! Errors:
  !!   fatalError if ID is not +ve
  !!
  subroutine setId(self, id)
    class(universe), intent(inout) :: self
    integer(shortInt), intent(in)  :: id
    character(100), parameter :: Here = 'setId (universe_inter.f90)'

    if (id <= 0) then
      call fatalError(Here, 'Id must be +ve. Is: '//numToChar(id))
    end if

    self % uniId = id

  end subroutine setId

  !!
  !! Get universe ID
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   universe ID.
  !!
  elemental function id(self)
    class(universe), intent(in) :: self
    integer(shortInt)           :: id

    id = self % uniID

  end function id

  !!
  !! Set universe index
  !!
  !! Args:
  !!   idx [in] -> Universe index (idx > 0)
  !!
  !! Errors:
  !!   fatalError if index is not +ve.
  !!
  subroutine setIdx(self, idx)
    class(universe), intent(inout) :: self
    integer(shortInt), intent(in)  :: idx
    character(100), parameter :: Here = 'setIdx (universe_inter.f90)'

    if (idx <= 0) then
      call fatalError(Here, 'Idx must be +ve. Is: '//numToChar(idx))
    end if

    self % uniIdx = idx

  end subroutine setIdx

  !!
  !! Read the data common to all universes from the dictionary
  !!
  !! Sets-up the base universe class with translations, ids and rotations (if any).
  !! Is used to avoid the code repeat in each subclass of the universe.
  !!
  !! Args:
  !!   dict [in] -> Dictionary with the universe definition
  !!
  subroutine setupBase(self, dict)
    class(universe), intent(inout)           :: self
    type(dictionary), intent(in)             :: dict
    integer(shortInt)                        :: id
    real(defReal), dimension(:), allocatable :: temp
    character(100), parameter :: Here = 'setupBase (universe_inter.f90)'

    ! Load basic data
    call dict % get(id, 'id')
    if (id <= 0) call fatalError(Here, 'Universe ID must be +ve. Is: ' // numToChar(id))
    call self % setId(id)

    ! Load origin
    if (dict % isPresent('origin')) then
      call dict % get(temp, 'origin')

      if (size(temp) /= 3) then
        call fatalError(Here, 'Origin must have size 3. Has: ' // numToChar(size(temp)))
      end if
      call self % setTransform(origin=temp)

    end if

    ! Load rotation
    if (dict % isPresent('rotation')) then
      call dict % get(temp, 'rotation')

      if (size(temp) /= 3) then
        call fatalError(Here, '3 rotation angles must be given. Has only: ' // numToChar(size(temp)))
      end if
      call self % setTransform(rotation=temp)
    end if

  end subroutine setupBase

  !!
  !! Set universe origin & rotation
  !!
  !! Note that rotation is defined by Euler angles with ZXZ convention
  !!
  !! Args:
  !!   origin [in]   -> Optional. 3D vector with the universe origin.
  !!   rotation [in] -> Optional. 3D vector with Euler ZXZ rotation angles [deg]. {phi, theta, psi}
  !!
  subroutine setTransform(self, origin, rotation)
    class(universe), intent(inout)                    :: self
    real(defReal), dimension(3), intent(in), optional :: origin
    real(defReal), dimension(3), intent(in), optional :: rotation

    if (present(origin)) then
      self % origin = origin

    end if

    if (present(rotation)) then
      if (.not.all(rotation == ZERO)) then ! Do not add matrix if there is no rotation
        self % rot = .true.
        call rotationMatrix(self % rotMat, rotation(1), rotation(2), rotation(3))
      end if
    end if

  end subroutine setTransform

  !!
  !! Enter from higher universe
  !!
  !! Sets:
  !!   - New position (r)
  !!   - New direction (dir)
  !!   - Rotation infor (isRotated + rotMat)
  !!   - Local cell (localID)
  !!   - Cell info (cellIdx)
  !!
  !! Does not set uniRootID !
  !! Applies all transformation when entering new universe with exception of cellOffset
  !! which needs to be applied to argument r outside this subroutine.
  !!
  !! Args:
  !!   new [out] -> New coordinates for the level.
  !!   r [in] -> Position affter cellOffset in upper univere is applied
  !!   u [in] -> Normalised direction (norm2(u) = 1.0)
  !!
  !! Note: Self is intent(inout), but if a state of the universe is to be changed
  !!   it is necessary to consider issues related to parallel calculations with shared
  !!   memory.
  !!
  subroutine enter(self, new, r, u)
    class(universe), intent(inout)          :: self
    type(coord), intent(out)                :: new
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u

    ! Set Info & Position
    new % r   = r
    new % dir = u
    new % uniIdx    = self % uniIdx
    new % isRotated = self % rot

    if (new % isRotated) then
      new % rotMat = self % rotMat
      new % r = matmul(self % rotMat, new % r)
      new % dir = matmul(self % rotMat, new % dir)
    end if

    ! Translate
    new % r = new % r - self % origin

    ! Find cell
    call self % findCell(new % localID, new % cellIdx, new % r, new % dir)

  end subroutine enter

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(universe), intent(inout) :: self

    self % uniIdx = 0
    self % uniId  = 0
    self % origin = ZERO
    self % rotMat = ZERO
    self % rot    = .false.

  end subroutine kill

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Utility Functions
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Convert a string to correct fill entry
  !!
  !! Material have names which are strings.
  !! Universes are identified by integer ID. Thus to refer to universe use
  !!   name = u<7>  ! For universe with ID 7
  !!   name = u<986> ! For universe with ID 986
  !!   ...
  !!
  !! Args:
  !!   name [in]  -> Character with the name of material or fill universe ID. Must not contain
  !!     any blanks!
  !!   mats [in]  -> Map of material names to matIdx
  !!   where [in] -> Name of caller procedure for better error message
  !!
  !! Result:
  !!   Filling. If material it is just matIdx. If universe it is -uniID (universe ID).
  !!
  function charToFill(name, mats, where) result(fill)
    character(nameLen), intent(in) :: name
    type(charMap), intent(in)      :: mats
    character(*), intent(in)       :: where
    integer(shortInt)              :: fill
    character(nameLen)             :: str
    integer(shortInt)              :: pos
    logical(defBool)               :: err
    integer(shortInt), parameter :: NOT_FOUND = -7

    ! Identify if mat or universe
    str = adjustl(name)
    if (str(1:2) == 'u<') then
      pos = index(str, '>', back=.true.)
    else
      pos = 0
    end if

    ! Convert to fill
    if (pos /= 0) then
      fill = charToInt(str(3:pos-1), error=err)
      if (err) call fatalError(where, 'Failed to convert '//trim(str)//' to universe ID')
      if (fill <= 0) call fatalError(where, 'Universe ID must be +ve is: '//numToChar(fill))
      fill = -fill

    else ! Convert to mat
      fill = mats % getOrDefault(name, NOT_FOUND)
      if (fill == NOT_FOUND) call fatalError(where, 'Unknown material: '//trim(name))
    end if

  end function charToFill

end module universe_inter
