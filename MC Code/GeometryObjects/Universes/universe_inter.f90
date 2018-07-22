module universe_inter

  use numPrecision
  use vector_class,  only : vector
  use coord_class,   only : coord
  use surface_inter, only : surfaceShelf
  use cell_class,    only : cellShelf

  implicit none
  private

  !!
  !! Universe is a collection of cells that spans the entire space
  !! It allows to find cellIdx and localId based on coordinates
  !! It determines distance to the next cell
  !! It performs cell-to-cell crossings thus:
  !!   -> special universes (pins, lattices) can use extra knowlage about likely neighbours
  !!
  type, public, abstract :: universe
    private
    integer(shortInt)          :: uniId
    real(defReal),dimension(3) :: offset = ZERO
  contains
    ! Build procedures
    procedure :: id
    procedure :: setId
    procedure :: setOffset

    ! Runtime procedures
    procedure                      :: enter
    procedure(findCell),deferred   :: findCell
    procedure(distance),deferred   :: distance
    procedure(cross),deferred      :: cross
    procedure(cellOffset),deferred :: cellOffset

  end type universe

  type,public,extends(universe) :: universeSlot
    private
    class(universe),allocatable :: slot
  contains
    ! Build procedures
    procedure :: id         => id_slot
    procedure :: setId      => setId_slot
    procedure :: setOffset  => setOffset_slot

    ! Run time procedures
    procedure :: enter      => enter_slot
    procedure :: findCell   => findCell_slot
    procedure :: distance   => distance_slot
    procedure :: cross      => cross_slot
    procedure :: cellOffset => cellOffset_slot

  end type universeSlot


  abstract interface
    !!
    !! Using the coordinates it finds a localID & cellIDx inside the universe
    !!
    subroutine findCell(self, coords, cShelf, sShelf)
      import :: universe   ,&
                coord      ,&
                cellShelf  ,&
                surfaceShelf
      class(universe), intent(in)    :: self
      type(coord),intent(inout)      :: coords
      type(cellShelf), intent(in)    :: cShelf
      type(surfaceShelf), intent(in) :: sShelf
    end subroutine findCell

    !!
    !! Returns distance to the next cell boundary in the universe
    !! Returns surfIdx of the surface beeing X-ed
    !! surfIdx can be < 0 (invalid)
    !! Invalid surfIdx indicate a surface that is private to the universe(not present on sShelf)
    !! Example of a private surface is face of a lattice cell
    !!
    subroutine distance(self, dist, surfIdx, coords ,cShelf, sShelf)
      import :: universe    ,&
                defReal     ,&
                shortInt    ,&
                coord       ,&
                cellShelf   ,&
                surfaceShelf
      class(universe), intent(in)    :: self
      real(defReal), intent(out)     :: dist
      integer(shortInt), intent(out) :: surfIdx
      type(coord), intent(in)        :: coords
      type(cellShelf), intent(in)    :: cShelf
      type(surfaceShelf),intent(in)  :: sShelf
    end subroutine distance


    !!
    !! Perform crossing inside the universe from current cell to next cell
    !! Assumes coords are at the surface being crossed (within surface_tol)
    !! IT DOES NOT PERFORM CHECK IF THE ABOVE ASSUMPTION IS VALID! (Be careful - MAK)
    !!
    subroutine cross(self, coords, surfIdx, cShelf, sShelf)
      import :: universe   ,&
                coord      ,&
                shortInt   ,&
                cellShelf  ,&
                surfaceShelf
      class(universe), intent(in)    :: self
      type(coord),intent(inout)      :: coords
      integer(shortInt), intent(in)  :: surfIdx
      type(cellShelf), intent(in)    :: cShelf
      type(surfaceShelf), intent(in) :: sShelf
    end subroutine cross

    !!
    !! Return offset for the current cell
    !! This is used when going into nested universe
    !! Total offset of the nested universe is :
    !!  cellOffset + nestedUniverse % offset
    !!
    !! Cell offset needs to be applied before envoing "enter" on the nested universe
    !!
    function cellOffset(self,coords) result(offset)
      import :: universe ,&
                coord    ,&
                defReal
      class(universe), intent(in)    :: self
      type(coord),intent(inout)      :: coords
      real(defReal),dimension(3)     :: offset
    end function cellOffset

  end interface

contains
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! universe procedures
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  !!
  !! Return id of the universe definition
  !!
  elemental function id(self)
    class(universe), intent(in) :: self
    integer(shortInt)           :: id

    id = self % uniID

  end function id

  !!
  !! Set universe Id
  !!
  elemental subroutine setId(self,id)
    class(universe), intent(inout) :: self
    integer(shortInt),intent(in)   :: id

    self % uniId = id

  end subroutine setId

  !!
  !! Set universe offset
  !!
  pure subroutine setOffset(self,offset)
    class(universe), intent(inout)        :: self
    real(defReal),dimension(3),intent(in) :: offset

    self % offset = offset

  end subroutine setOffset
    
  !!
  !! Enter universe
  !! Apply co-ordinate transformation and find cell
  !! Set localID and cellIdx in coords
  !!
  subroutine enter(self, coords, cShelf, sShelf)
    class(universe), intent(in)    :: self
    type(coord),intent(inout)      :: coords
    type(cellShelf), intent(in)    :: cShelf
    type(surfaceShelf), intent(in) :: sShelf

    ! Transform coordinates
    coords % r = coords % r - self % offset

    ! Find cell
    call self % findCell(coords, cShelf, sShelf)

  end subroutine enter

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! universeSlot procedures
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Get id of universe in the slot
  !!
  elemental function id_slot(self) result(id)
    class(universeSlot), intent(in) :: self
    integer(shortInt)               :: id

    id = self % slot % id()

  end function id_slot

  !!
  !! Set id of universe in the slot
  !!
  elemental subroutine setId_slot(self,id)
    class(universeSlot), intent(inout) :: self
    integer(shortInt),intent(in)       :: id

    call self % slot % setId(id)
    self % uniId = id

  end subroutine setId_slot

  !!
  !! Set offset of universe in the slot
  !!
  pure subroutine setOffset_slot(self,offset)
    class(universeSlot), intent(inout)    :: self
    real(defReal),dimension(3),intent(in) :: offset

    call self % slot % setOffset(offset)
    self % offset = offset

  end subroutine setOffset_slot

  !!
  !! Enter universe in the slot
  !!
  subroutine enter_slot(self, coords, cShelf, sShelf)
    class(universeSlot), intent(in)  :: self
    type(coord),intent(inout)        :: coords
    type(cellShelf), intent(in)      :: cShelf
    type(surfaceShelf), intent(in)   :: sShelf

    call self % slot % enter(coords, cShelf, sShelf)

  end subroutine enter_slot

  !!
  !! Find cell in the universe in the slot
  !!
  subroutine findCell_slot(self, coords, cShelf, sShelf)
    class(universeSlot), intent(in)  :: self
    type(coord),intent(inout)        :: coords
    type(cellShelf), intent(in)      :: cShelf
    type(surfaceShelf), intent(in)   :: sShelf

    call self % slot % findCell(coords, cShelf, sShelf)

  end subroutine findCell_slot

  !!
  !! Find distance to next cell in the universe in the slot
  !!
  subroutine distance_slot(self, dist, surfIdx, coords, cShelf, sShelf)
    class(universeSlot), intent(in)  :: self
    real(defReal), intent(out)       :: dist
    integer(shortInt), intent(out)   :: surfIdx
    type(coord), intent(in)          :: coords
    type(cellShelf), intent(in)      :: cShelf
    type(surfaceShelf),intent(in)    :: sShelf

    call self % slot% distance(dist, surfIdx, coords, cShelf, sShelf)

  end subroutine distance_slot

  !!
  !! Cross surface in the universe in the slot
  !!
  subroutine cross_slot(self, coords, surfIdx, cShelf, sShelf)
    class(universeSlot), intent(in)  :: self
    type(coord),intent(inout)        :: coords
    integer(shortInt), intent(in)    :: surfIdx
    type(cellShelf), intent(in)      :: cShelf
    type(surfaceShelf), intent(in)   :: sShelf

    call self % slot % cross(coords, surfIdx, cShelf, sShelf)

  end subroutine cross_slot

  !!
  !! Give offset for current cell in the universe in the slot
  !!
  function cellOffset_slot(self,coords) result(offset)
    class(universeSlot), intent(in) :: self
    type(coord),intent(inout)       :: coords
    real(defReal),dimension(3)      :: offset

    offset = self % slot % cellOffset(coords)

  end function cellOffset_slot

end module universe_inter
