module cell_class
  use numPrecision
  use universalVariables
  use genericProcedures, only : fatalError, linFind, targetNotFound
  use vector_class,      only : vector
  use surface_inter,     only  : surface, surfaceShelf

  implicit none
  private

  !! Local Parameters
  integer(shortInt),parameter    :: DEF_STRIDE = 20

  !!
  !! Geometry Cell
  !! Cell is an intersection of surface halfspaces
  !! +ve and -ve halfspaces are indentified with a sign of surfIdx
  !!
  !! NOTE: Definition of a cell tells us nothing about its content. It only defines region
  !! of space the cell occupies. Content of the cell is obtained from fill Array by providing
  !! uniqueID of a cell
  !!
  type, public :: cell
    private
    integer(shortInt)                            :: cellId      ! unique user-defined ID
    integer(shortInt)                            :: cellIdx     ! Index of the cell on cellShelf
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

    ! Runtime procedure
    procedure          :: isInside
    procedure          :: distance

  end type cell


  type, public :: cellShelf
    private
    type(cell),dimension(:),allocatable, public :: shelf
    integer(shortInt)                           :: Nmax   = 0
    integer(shortInt)                           :: N      = 0
    integer(shortInt)                           :: stride = DEF_STRIDE
    integer(shortInt)                           :: nextID = 1
  contains
    procedure :: init       => init_shelf
    procedure :: kill       => kill_shelf
    procedure :: trimSize   => trimSize_shelf

    procedure :: addUnique  => addUnique_shelf
    !procedure :: getOrAdd   => getOrAdd_shelf
    procedure :: getIdx     => getIdx_shelf


    ! Private procedures
    procedure, private :: grow    => grow_shelf
    procedure, private :: resize  => resize_shelf
    procedure, private :: freeID  => freeId_shelf
  end type
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
  subroutine init(self, surfIds, id, cellIdx, surfShelf)
    class(cell), intent(inout)                   :: self
    integer(shortInt),dimension(:),intent(in)    :: surfIDs
    integer(shortInt), intent(in)                :: id
    integer(shortInt), intent(in)                :: cellIdx
    type(surfaceShelf), intent(in)               :: surfShelf
    character(100), parameter :: Here ='init (cell_class.f90)'

    ! Load surface IDXs
    self % surfaces = surfShelf % getIdx( abs(surfIDs) )

    ! Verify that all IDs were present
    if(any( self % surfaces == targetNotFound)) then
      print *, surfIds
      print *, self % surfaces
      call fatalError(Here,'One or more of the surface IDs was not found on the shelf')
    end if

    ! Correct sign of surfaces
    self % surfaces = sign(self % surfaces, surfIDs)

    ! Load cellIdx and Id
    self % cellId  = id
    self % cellIdx = cellIdx

  end subroutine init

  !!
  !! Move cell from RHS to LHS
  !! Moves allocation of "surface" component
  !!
  elemental subroutine moveAllocFrom(LHS,RHS)
    class(cell), intent(inout) :: LHS
    type(cell), intent(inout)  :: RHS

    ! Copy components
    LHS % cellId  = RHS % cellId
    LHS % cellIdx = RHS % cellIdx

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
    self % cellIdx = 0
    if(allocated(self % surfaces)) deallocate(self % surfaces)

  end subroutine kill

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
        idx  = testIdx

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
  !! Adds allocatable surface to the shelf
  !! Surf needs to be allocated upon entry
  !! Deallocates surface in the process
  !! Throws error if surface with the same definition or id is already present
  !! Returns surfIdx in the surfaceShelf
  !!
  subroutine addUnique_shelf(self,newCell,cellIdx)
    class(cellShelf), intent(inout)          :: self
    type(cell),intent(inout)                 :: newCell
    integer(shortInt), intent(out)           :: cellIdx
    integer(shortInt)                        :: N
    character(100),parameter :: Here = 'addUnique_shelf (cell_class.f90)'

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

        else if ( any( storedCells .same. newCell )) then          ! Definition is not unique
          call fatalError(Here,'Cell with the same definition is already present')

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

  end subroutine addUnique_shelf

  !!
  !! Returns surfId and surfIdx of the surface with the same definition
  !! Surf needs to be allocated upon entry
  !! If there is no surface with the same definition adds it the shelf
  !! It completly ignores id of surf dummy argument
  !! surf if deallocated upon exit
  !!
  subroutine getOrAdd_shelf(self,newCell,cellId,cellIdx)
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
    call newCell % setId(surfId)
    call self % addUnique(surf,surfIdx)

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
