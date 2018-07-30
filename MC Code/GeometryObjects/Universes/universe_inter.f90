module universe_inter

  use numPrecision
  use genericProcedures,  only : fatalError, linFind
  use vector_class,       only : vector
  use coord_class,        only : coord
  use surface_inter,      only : surfaceShelf
  use cell_class,         only : cellShelf

  implicit none
  private

  !! Local Parameters
  integer(shortInt),parameter    :: DEF_STRIDE = 20

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
    procedure                      :: id
    procedure                      :: setId
    procedure                      :: setOffset
    procedure(cellIdx), deferred   :: cellIdx

    ! Runtime procedures
    procedure                      :: enter
    procedure(findCell),deferred   :: findCell
    procedure(distance),deferred   :: distance
    procedure(cross),deferred      :: cross
    procedure(cellOffset),deferred :: cellOffset

  end type universe

  !!
  !! Slot to store diffrent polymorphic universes in the single array
  !!
  type,public,extends(universe) :: universeSlot
    private
    class(universe),allocatable :: slot
  contains
    ! Loading into slot
    generic           :: assignment(=) => copy_slot
    procedure         :: moveAllocFrom => moveAllocFrom_slot
    procedure         :: load          => load_slot
    procedure,private :: copy_slot

    ! Build procedures
    procedure :: id         => id_slot
    procedure :: setId      => setId_slot
    procedure :: setOffset  => setOffset_slot
    procedure :: cellIdx    => cellIDx_slot

    ! Run time procedures
    procedure :: enter      => enter_slot
    procedure :: findCell   => findCell_slot
    procedure :: distance   => distance_slot
    procedure :: cross      => cross_slot
    procedure :: cellOffset => cellOffset_slot

  end type universeSlot

  !!
  !! Stores number of universes
  !! Public interface:
  !!  shelf     -> array of universes
  !!  add       -> adds a universe with unique id to Shelf
  !!  getOrAdd  -> returns ID & Idx of provided universe, If it isn't present it adds it to the shelf
  !!  getIdx    -> returns Idx of universe with ID. returns "targetNotFound" if ID is not present
  !!  trimSize  -> Resizes array to contain no empty entries
  !!  init      -> initialises with "maximum size" and "stride"
  !!  kill      -> returns to uninitialised state
  !!
  type, public :: universeShelf
    private
    type(universeSlot),dimension(:),allocatable, public :: shelf
    integer(shortInt)                                   :: Nmax   = 0
    integer(shortInt)                                   :: N      = 0
    integer(shortInt)                                   :: stride = DEF_STRIDE
    integer(shortInt)                                   :: nextID = 1
  contains
    procedure :: init       => init_shelf
    procedure :: kill       => kill_shelf
    procedure :: trimSize   => trimSize_shelf
    procedure :: size       => size_shelf

    procedure :: add        => add_shelf
    procedure :: getIdx     => getIdx_shelf

    ! Private procedures
    procedure, private :: grow    => grow_shelf
    procedure, private :: resize  => resize_shelf

  end type universeShelf


  abstract interface
    !!
    !! Given cell localID returtns cellIdx
    !!
    function cellIdx(self, localID)
      import :: universe, &
                shortInt
      class(universe), intent(in)   :: self
      integer(shortInt), intent(in) :: localID
      integer(shortInt)             :: cellIdx
    end function cellIdx

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
  !! Load allocatable surface on RHS into LHS slot
  !! Be carefull about loading slots into slots.
  !! It will work but the chain will hurt performance
  !!
  subroutine load_slot(LHS,RHS)
    class(universeSlot), intent(inout)         :: LHS
    class(universe),allocatable, intent(inout) :: RHS

    if(allocated(LHS % slot)) deallocate (LHS % slot)
    call move_alloc(RHS , LHS % slot)

    LHS % uniID  = LHS % slot % id()
    LHS % offset = LHS % slot % offset

  end subroutine load_slot

  !!
  !! Copy RHS into LHS
  !! Be carefull about loading slots into slots.
  !! It will work but the chain will hurt performance
  !!
  elemental subroutine copy_slot(LHS,RHS)
    class(universeSlot),intent(inout) :: LHS
    class(universe), intent(in)       :: RHS

    if(allocated(LHS % slot)) deallocate (LHS % slot)
    allocate(LHS % slot, source = RHS)

    LHS % uniID  = LHS % slot % id()
    LHS % offset = LHS % slot % offset

  end subroutine copy_slot

  !!
  !! Move allocation from RHS slot into LHS slot
  !!
  elemental subroutine moveAllocFrom_slot(LHS,RHS)
    class(universeSlot), intent(inout) :: LHS
    class(universeSlot), intent(inout) :: RHS

    if(allocated(LHS % slot)) deallocate (LHS % slot)
    call move_alloc(RHS % slot, LHS % slot)

    LHS % uniID  = LHS % slot % id()
    LHS % offset = LHS % slot % offset

  end subroutine moveAllocFrom_slot

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
  !! Return cellIdx given local cellId
  !!
  function cellIdx_slot(self, localID)
    class(universeSlot), intent(in) :: self
    integer(shortInt), intent(in)   :: localID
    integer(shortInt)               :: cellIdx_slot

    cellIdx_slot = self % slot % cellIdx(localID)

  end function cellIdx_slot

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

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! universeShelf procedures
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  !!
  !! Initialise shelf
  !! Give initial size and optionaly stride to override default
  !!
  subroutine init_shelf(self,N,stride)
    class(universeShelf), intent(inout)    :: self
    integer(shortInt), intent(in)          :: N
    integer(shortInt), optional,intent(in) :: stride

    ! Return to uninitialised state
    call self % kill()

    ! Allocate storage space
    allocate( self % shelf(N))

    ! Assign constants
    self % Nmax = N
    if(present(stride)) self % stride = stride

  end subroutine init_shelf

  !!
  !! Deallocate storage space and returns shelf to uninitialised state
  !!
  subroutine kill_shelf(self)
    class(universeShelf), intent(inout) :: self

    if(allocated(self % shelf)) deallocate(self % shelf)
    self % Nmax   = 0
    self % N      = 0
    self % stride = DEF_STRIDE
    self % nextID = 1

  end subroutine kill_shelf

  !!
  !! Set shelf size to the number of stored entries
  !!
  subroutine trimSize_shelf(self)
    class(universeShelf), intent(inout) :: self

    call self % resize(self % N)

  end subroutine trimSize_shelf

  !!
  !! Returns current size of shelf filled with non-empty elements
  !!
  function size_shelf(self) result(N)
    class(universeShelf), intent(in) :: self
    integer(shortInt)                :: N

    N = self % N

  end function size_shelf

  !!
  !! Adds an universe to the shelf
  !! Universe needs to be initialised
  !! Kills universe in the process
  !! Throws error if univerese with the same id is already present
  !! Returns uniIdx in the universeShelf
  !!
  subroutine add_shelf(self,newUni,uniIdx)
    class(universeShelf), intent(inout)         :: self
    class(universe), allocatable ,intent(inout) :: newUni
    integer(shortInt), intent(out)              :: uniIdx
    integer(shortInt)                           :: N
    character(100),parameter :: Here = 'add_shelf (universe_inter.f90)'

    if(.not.allocated(newUni )) then
      call fatalError(Here, 'Provided universe to store was not initialised')
    end if

    ! Load surface counter for convenience
    N = self % N

    if( N == 0) then ! First surface is  being loaded
      call self % shelf(1) % load(newUni)
      uniIdx = 1

    else
      ! There are already surfaces present (need to check uniqueness)
      associate ( storedUnis => self % shelf(1:N) )

        if ( any( storedUnis % id() == newUni % id() )) then        ! ID is already present
          print *, storedUnis % id()
          print *, newUni % id()
          call fatalError(Here,'Id is already present in the universe shelf')

        end if
      end associate

      ! Load surface
      call self % shelf(N+1) % load(newUni)
      uniIdx = N +1

    end if

    ! Increment surface counter
    self % N = N + 1

    ! Grow storage if needed
    call self % grow()

  end subroutine add_shelf

  !!
  !! Return universe idx given universe ID
  !! Returns "targetNotFound" if id is not present
  !!
  elemental function getIdx_shelf(self,id) result(idx)
    class(universeShelf), intent(in) :: self
    integer(shortInt), intent(in)    :: id
    integer(shortInt)                :: idx
    integer(shortInt)                :: i

    associate( storedUniss => self % shelf(1:self % N))
      idx = linFind(storedUniss % id(), id)

    end associate
  end function getIdx_shelf

  !!
  !! Increases the size of the shelf by stride if N >= Nmax
  !!
  subroutine grow_shelf(self)
    class(universeShelf), intent(inout) :: self
    integer(shortInt)                   :: Nmax

    ! Return if shelf does not need to grow
    if ( self % N < self % Nmax) then
      return
    end if

    ! Calculate new maximum size
    Nmax = self % Nmax + self % stride
    call self % resize(Nmax)

  end subroutine grow_shelf

  !!
  !! Set size of the shelf to Nmax >= self % N
  !!
  subroutine resize_shelf(self,Nmax)
    class(universeShelf), intent(inout)          :: self
    integer(shortInt), intent(in)                :: Nmax
    type(universeSlot),dimension(:),allocatable  :: tempShelf
    character(100), parameter :: Here ='resize (cell_class.f90)'

    if( Nmax < self % N) call fatalError(Here,'Resizeing to size smaller then number of entries')

    ! Allocate new storage space
    allocate(tempShelf(Nmax) )

    ! Copy Old enteries without reallocation
    associate ( new => tempShelf(1:self % N), &
                old => self % shelf(1:self % N) )

      call new % moveAllocFrom( old)
    end associate

    ! Replace shelf and update Nmax
    call move_alloc(tempShelf, self % shelf)
    self % Nmax = Nmax

  end subroutine resize_shelf

end module universe_inter
