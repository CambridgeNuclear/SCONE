module cell_class
  use numPrecision
  use universalVariables
  use genericProcedures, only : fatalError, linFind, targetNotFound, numToChar, hasDuplicates
  use vector_class,      only : vector
  use dictionary_class,  only : dictionary
  use surface_inter,     only : surface, surfaceSlot, surfaceShelf

  implicit none
  private

  !! Local Parameters
  integer(shortInt),parameter    :: DEF_STRIDE = 20

  !!
  !! Public function to get cell from dict and surface shelf
  !!
  public :: new_cell

  !!
  !! Geometry Cell
  !! Cell is an intersection of surface halfspaces
  !! +ve and -ve halfspaces are indentified with a sign of surfIdx
  !!
  !! NOTE:
  !!   1) Definition of a cell tells us nothing about its content. It only defines region
  !!   of space the cell occupies. Content of the cell is obtained from fill Array by providing
  !!   uniqueID of a cell
  !!   2) Cells stores surfaces in order given by the user. This makes some checks for the same
  !!   definition mopre convoluted (permutations of surfaces need to be detected), but allows
  !!   a bit of hand optimisation when checking for halfspace (start with a surface particle is
  !!   least likley to be).
  !!
  type, public :: cell
    private
    integer(shortInt)                            :: cellId      ! unique user-defined ID
    integer(shortInt),dimension(:),allocatable   :: surfaces    ! +/- surfIDXs for +/- halfspace

  contains
    ! Operators and assignments
    generic   :: operator(.same.) => same_cell

    ! Build procedures
    procedure          :: init
    procedure          :: moveAllocFrom
    procedure          :: id
    procedure, private :: same_cell
    procedure          :: kill
    procedure          :: getDef

    ! Runtime procedure
    procedure          :: isInside
    procedure          :: distance

  end type cell

  !!
  !! Stores number of cells
  !! Public interface:
  !!  shelf     -> array of cells
  !!  add       -> adds a cell with unique id to Shelf
  !!  getOrAdd  -> returns ID & Idx of provided cell, If it isn't present it adds it to the shelf
  !!  getIdx    -> returns Idx of cell with ID. returns "targetNotFound" if ID is not present
  !!  trimSize  -> Resizes array to contain no empty entries
  !!  init      -> initialises with "maximum size" and "stride"
  !!  kill      -> returns to uninitialised state
  !!
  type, public :: cellShelf
    !private
    type(cell),dimension(:),allocatable, public :: shelf
    integer(shortInt)                           :: Nmax   = 0
    integer(shortInt)                           :: N      = 0
    integer(shortInt)                           :: stride = DEF_STRIDE
    integer(shortInt)                           :: nextID = 1
  contains
    procedure :: init       => init_shelf
    procedure :: kill       => kill_shelf
    procedure :: trimSize   => trimSize_shelf
    procedure :: size       => size_shelf

    procedure :: add        => add_shelf
    procedure :: getOrAdd   => getOrAdd_shelf
    procedure :: getIdx     => getIdx_shelf

    ! Private procedures
    procedure, private :: grow    => grow_shelf
    procedure, private :: resize  => resize_shelf
    procedure, private :: freeID  => freeId_shelf

  end type cellShelf

contains
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! cell procedures
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  !!
  !! Initialise the cell
  !! Provide vector of surface IDs. Let sign determine halfspace
  !! Give id and idx of the cell on the callShelf
  !! Give surfShelf to translate surfId to surfIdx
  !!
  subroutine init(self, surfIds, id, surfShelf)
    class(cell), intent(inout)                   :: self
    integer(shortInt),dimension(:),intent(in)    :: surfIDs
    integer(shortInt), intent(in)                :: id
    type(surfaceShelf), intent(in)               :: surfShelf
    character(100), parameter :: Here ='init (cell_class.f90)'

    ! Load surface IDXs
    self % surfaces = surfShelf % getIdx( abs(surfIDs) )

    ! Verify that all IDs were present
    if(any( self % surfaces == targetNotFound)) then
      print *, "Surface IDs :",surfIds
      print *, "Surface IDXs:",self % surfaces
      call fatalError(Here,'One or more of the surface IDs was not found on the shelf in cell: ' &
                            // numToChar(id) )
    end if

    ! Correct sign of surfaces
    self % surfaces = sign(self % surfaces, surfIDs)

    ! Load cellId
    self % cellId  = id

    ! Check that no surfaces are repeated in surfIds
    if(hasDuplicates(abs(self % surfaces))) then
      call fatalError(Here,'There are repeated surfaces in definition of cell: '//numToChar(id))
    end if

  end subroutine init

  !!
  !! Return a cell from a dictionary and surface shelf
  !!
  function new_cell(dict,sShelf) result(new)
    class(dictionary), intent(in)              :: dict
    type(surfaceShelf), intent(inout)          :: sShelf
    type(cell)                                 :: new
    integer(shortInt)                          :: id
    integer(shortInt),dimension(:),allocatable :: surfIds

    ! Obtain input
    call dict % get(id,'id')
    call dict % get(surfIds,'surfaces')

    ! Initialise
    call new % init(surfIds, id, sShelf)

  end function new_cell

  !!
  !! Move cell from RHS to LHS
  !! Moves allocation of "surface" component
  !!
  elemental subroutine moveAllocFrom(LHS,RHS)
    class(cell), intent(inout) :: LHS
    type(cell), intent(inout)  :: RHS

    ! Copy components
    LHS % cellId  = RHS % cellId

    ! Move allocation
    call move_alloc(RHS % surfaces, LHS % surfaces)

  end subroutine moveAllocFrom

  !!
  !! Return ID of the cell
  !!
  elemental function id(self)
    class(cell), intent(in) :: self
    integer(shortInt)       :: id

    id = self % cellId

  end function id

  !!
  !! Compares definition of two cells and returns .true. if they are the same
  !!
  elemental function same_cell(LHS,RHS) result(same)
    class(cell),intent(in) :: LHS
    class(cell),intent(in) :: RHS
    logical(defBool)       :: same
    integer(shortInt)      :: i, testIdx

    ! Compare size of surface component
    same = size(LHS % surfaces) == size(RHS % surfaces)

    ! Compare individual entries if size is the same
    ! *** Need to account for possible permutations in definitions
    if(same) then
      do i=1,size(LHS % surfaces)
        testIdx = linFind(RHS % surfaces, LHS % surfaces(i))

        if (testIdx == targetNotFound) then
          same = .false.
          return
        end if

      end do
    end if

  end function same_cell

  !!
  !! Returns cell to uninitialised state
  !!
  elemental subroutine kill(self)
    class(cell), intent(inout) :: self

    self % cellId  = 0
    if(allocated(self % surfaces)) deallocate(self % surfaces)

  end subroutine kill

  !!
  !! Returns array of surfaces defining cell with +/-ve signs for halfspaces
  !!
  function getDef(self) result(def)
    class(cell), intent(in)                    :: self
    integer(shortInt),dimension(:),allocatable :: def

    def = self % surfaces

  end function getDef

  !!
  !! Checks whether a point occupies a cell by examining each surface in turn
  !! Returns TRUE if point is in cell, FALSE if point is not
  !!
  elemental function isInside(self, r, u, surfShelf) result(isIt)
    class(cell), intent(in)                 :: self
    type(vector), intent(in)                :: r
    type(vector), intent(in)                :: u
    type(surfaceShelf), intent(in)          :: surfShelf
    logical(defBool)                        :: isIt
    integer(shortInt)                       :: surfIdx, i
    logical(defBool)                        :: halfspace, sense

    ! Need only check that halfspaces are satisfied: if not, the point is outside the cell
    ! Check each cell surface to ensure point is within its halfspace
     do i=1,size(self % surfaces)
       surfIdx   = abs(self % surfaces(i))

       ! Check that halfspace of the surface is equivalent to one in definition
       halfspace = surfShelf % shelf(surfIdx) % halfspace(r,u)
       sense     = self % surfaces(i) > 0

       if ( halfspace .neqv. sense) then
         isIt = .false.
         return

       end if
     end do

     ! If reached here it means point is inside cell
     isIt = .true.

  end function isInside

  !!
  !! Find the shortest positive distance to the boundary of the cell
  !!
  elemental subroutine distance(self, dist, idx, r, u, surfShelf)
    class(cell), intent(in)          :: self
    real(defReal), intent(out)       :: dist
    integer(shortInt), intent(out)   :: idx
    type(vector), intent(in)         :: r
    type(vector), intent(in)         :: u
    type(surfaceShelf), intent(in)   :: surfShelf
    real(defReal)                    :: testDist
    integer(shortInt)                :: i, testIdx, surfIdx

    dist = INFINITY
    idx  = 0

    ! Search through all surfaces to find the minimum distance to the boundary
    ! Should not have to check for negative distances: these are set to INFINITY
    ! in the distance routines
    do i = 1, size(self % surfaces)

      ! Obtain test distance and corresponding surface index
      surfIdx = abs(self % surfaces(i))
      call surfShelf % shelf(surfIdx) % distance(testDist, testIdx, r, u)

      ! Retain smallest distance & corresponding surfaceIdx
      ! *** We use testIdx becouse compond surfaces can return index of other surface
      if (testDist < dist) then
        dist = testDist
        idx  = surfIdx
        ! For now return index of the global surface
        !idx  = testIdx

      end if
    end do
  end subroutine distance

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! cellShelf procedures
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  !! Initialise shelf
  !! Give initial size and optionaly stride to override default
  !!
  subroutine init_shelf(self,N,stride)
    class(cellShelf), intent(inout)        :: self
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
    class(cellShelf), intent(inout) :: self

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
    class(cellShelf), intent(inout) :: self

    call self % resize(self % N)

  end subroutine trimSize_shelf

  !!
  !! Returns current size of shelf filled with non-empty elements
  !!
  function size_shelf(self) result(N)
    class(cellShelf), intent(in) :: self
    integer(shortInt)            :: N

    N = self % N

  end function size_shelf

  !!
  !! Adds acell to the shelf
  !! Cell needs to be initialised
  !! Kills cell in the process
  !! Throws error if cell with the same id is already present
  !! Returns cellIdx in the cellShelf
  !!
  subroutine add_shelf(self,newCell,cellIdx)
    class(cellShelf), intent(inout)          :: self
    type(cell),intent(inout)                 :: newCell
    integer(shortInt), intent(out)           :: cellIdx
    integer(shortInt)                        :: N
    character(100),parameter :: Here = 'add_shelf (cell_class.f90)'

    if(.not.allocated(newCell % surfaces )) then
      call fatalError(Here, 'Provided cell to store was not initialised')
    end if

    ! Load surface counter for convenience
    N = self % N

    if( N == 0) then ! First surface is  being loaded
      call self % shelf(1) % moveAllocFrom(newCell)
      cellIdx = 1

    else
      ! There are already surfaces present (need to check uniqueness)
      associate ( storedCells => self % shelf(1:N) )

        if ( any( storedCells % id() == newCell % id() )) then        ! ID is already present
          call fatalError(Here,'Id is already present in the cell shelf')

        end if
      end associate

      ! Load surface
      call self % shelf(N+1) % moveAllocFrom(newCell)
      cellIdx = N +1

    end if

    ! Increment surface counter
    self % N = N + 1

    ! Grow storage if needed
    call self % grow()

  end subroutine add_shelf

  !!
  !! Returns cellId and cellIdx of the cell with the same definition
  !! Cell needs to be initialised on entry
  !! If there is no cell with the same definition adds it the shelf
  !! It completly ignores id of cell dummy argument
  !! cell if uninitialised upon exit
  !!
  subroutine getOrAdd_shelf(self, newCell, cellId, cellIdx)
    class(cellShelf), intent(inout)        :: self
    type(cell), intent(inout)              :: newCell
    integer(shortInt), intent(out)         :: cellId
    integer(shortInt), intent(out)         :: cellIdx
    integer(shortInt)                      :: i
    character(100),parameter :: Here = 'getOrAdd_shelf (cell_class.f90)'

    if(.not.allocated(newCell % surfaces )) then
      call fatalError(Here, 'Provided cell to store was not initialised')
    end if

    ! Find surface with the same definition
    do i=1,self % N
      if(newCell .same. self % shelf(i)) then ! Return Idx & ID of existing cell
        cellId  = self % shelf(i) % id()
        cellIdx = i
        call newCell % kill()
        return

      end if
    end do

    ! Surface was not found set a next free Id and add to shelf
    cellId = self % freeId()
    newCell % cellId = cellId
    call self % add(newCell, cellIdx)

  end subroutine getOrAdd_shelf

  !!
  !! Return cell idx given cell ID
  !! Returns "targetNotFound" if id is not present
  !!
  elemental function getIdx_shelf(self,id) result(idx)
    class(cellShelf), intent(in)  :: self
    integer(shortInt), intent(in) :: id
    integer(shortInt)             :: idx
    integer(shortInt)             :: i

    associate( storedCells => self % shelf(1:self % N))
      idx = linFind(storedCells % id(), id)

    end associate
  end function getIdx_shelf


  !!
  !! Increases the size of the shelf by stride if N >= Nmax
  !!
  subroutine grow_shelf(self)
    class(cellShelf), intent(inout)            :: self
    integer(shortInt)                          :: Nmax

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
    class(cellShelf), intent(inout)      :: self
    integer(shortInt), intent(in)        :: Nmax
    type(cell),dimension(:),allocatable  :: tempShelf
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

  !!
  !! Returns smallest +ve free id in the shelf
  !!
  function freeId_shelf(self) result(id)
    class(cellShelf), intent(inout) :: self
    integer(shortInt)               :: id

    associate( storedCells => self % shelf(1:self % N) )

      do
        id = self % nextId

        ! Check if the proposed id is unique
        if ( all(storedCells % id() /= id )) return

        ! Increment test id if current nextId is occupied
        self % nextId = self % nextId + 1

      end do
    end associate

  end function freeId_shelf

end module cell_class
