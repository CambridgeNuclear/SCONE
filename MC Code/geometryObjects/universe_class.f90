!
! The universe class contains collections of cells
! Universes extend infinitely
!
! The universe can be searched to find a cell.
! The universe also provides co-ordinate transformations relative to its parent universe
!

module universe_class
  use numPrecision
  use genericProcedures
  use surface_class
  use cell_class

  implicit none
  private

  type, public :: universe
    type(cell_ptr), dimension(:), allocatable :: cells     ! cells contained in universe
    integer(shortInt) :: numCells                          ! the number of constituent cells
    real(defReal), dimension(3) :: offset                  ! x,y,z offset from parent universe
    integer(shortInt) :: id                                ! unique universe ID
    integer(shortInt) :: geometryInd
    logical(defBool) :: rootUni = .false.                  ! Is the universe the root universe?
    character(100) :: name = ""
  contains
    procedure :: init                                      ! initialise universe
    procedure :: whichCell                                 ! Identify which cell a particle occupies
  end type universe

  type, public :: universe_ptr
    class(universe), pointer :: ptr
  contains
    procedure :: init => init_ptr
    procedure :: whichCell => whichCell_ptr
    procedure :: cells => cells_ptr
    procedure :: offset => offset_ptr
    procedure :: name => name_ptr
    procedure :: geometryInd => geometryInd_ptr
    procedure :: numCells => numCells_ptr
    procedure :: associated => associated_ptr
    procedure :: kill
    generic   :: assignment(=) => universe_ptr_assignment, universe_ptr_assignment_target
    procedure,private :: universe_ptr_assignment
    procedure,private :: universe_ptr_assignment_target
  end type universe_ptr

contains

  !!
  !! Initialise the universe
  !!
  subroutine init(self, offset, cells, id, geometryInd, isRoot, name)
    class(universe), intent(inout) :: self
    real(defReal), dimension(3), intent(in) :: offset
    class(cell_ptr), dimension(:), intent(in) :: cells
    logical(defBool), optional, intent(in) :: isRoot
    integer(shortInt), intent(in) :: id
    integer(shortInt), intent(in) :: geometryInd
    character(*), optional, intent(in) :: name

    self % numCells = size(cells)
    allocate(self % cells(self % numCells))
    self % offset = offset
    self % cells = cells
    self % id = id
    self % geometryInd = geometryInd
    if(present(isRoot)) then
      self % rootUni = isRoot
    else
      self % rootUni = .false.
    end if
    if(present(name)) self % name = name
  end subroutine init

  !!
  !! Search cells to find which cell a particle occupies
  !! Returns the index of the cell in the cells array
  !!
  function whichCell(self, r, u) result(c)
    class(universe), intent(in) :: self
    real(defReal), dimension(3), intent(in) :: r, u
    type(cell_ptr) :: c
    integer(shortInt) :: i
    logical(defBool) :: isInside

    ! Loop through each cell within the universe, terminate the search if found
    do i=1, self % numCells
      isInside = self % cells(i) % insideCell(r, u)
      if (isInside) then
        c = self % cells(i)
        return
      end if
    end do

    call fatalError('whichCell, universe','Point could not be found')
    ! Dummy statement to avoid compiler warning. Will(Should) never be executed...
    c = self % cells(1)
  end function whichCell

!
! Pointer wrapper functions
!
!

  !!
  !! Initialise the root universe
  !!
  subroutine init_ptr(self, offset, cells, id, geometryInd, isRoot, name)
    class(universe_ptr), intent(inout) :: self
    real(defReal), dimension(3), intent(in) :: offset
    class(cell_ptr), dimension(:), intent(in) :: cells
    integer(shortInt), intent(in) :: id
    integer(shortInt), intent(in) :: geometryInd
    logical(defBool), optional, intent(in) :: isRoot
    character(*), optional, intent(in) :: name
    call self % ptr % init(offset, cells, id, geometryInd, isRoot, name)
  end subroutine init_ptr


  !!
  !! Search cells to find which cell a particle occupies
  !! Returns the index of the cell in the cells array
  !!
  function whichCell_ptr(self, r, u) result(c)
    class(universe_ptr), intent(in) :: self
    real(defReal), dimension(3), intent(in) :: r, u
    type(cell_ptr) :: c
    c = self % ptr % whichCell(r,u)
  end function whichCell_ptr

  !!
  !! Returns a pointer to a cell within the universe
  !!
  function cells_ptr(self,i)result(c)
    class(universe_ptr), intent(in) :: self
    integer(shortInt), intent(in) :: i
    type(cell_ptr) :: c
    c = self % ptr % cells(i)
  end function cells_ptr

  !!
  !! Returns the offset of the universe
  !!
  function offset_ptr(self)result(offset)
    class(universe_ptr), intent(in) :: self
    real(defReal),dimension(3) :: offset
    offset = self % ptr % offset
  end function offset_ptr

  !!
  !! Returns the name of the universe pointed to
  !!
  function name_ptr(self)result(name)
    class(universe_ptr), intent(in) :: self
    character(100) :: name
    name = self % ptr % name
  end function name_ptr

  !!
  !! Returns the geometry index of the universe pointed to
  !!
  function geometryInd_ptr(self)result(ind)
    class(universe_ptr), intent(in) :: self
    integer(shortInt) :: ind
    ind = self % ptr % geometryInd
  end function geometryInd_ptr

  !!
  !! Returns the number of cells within the universe
  !!
  function numCells_ptr(self)result(num)
    class(universe_ptr), intent(in) :: self
    integer(shortInt) :: num
    num = self % ptr % numCells
  end function numCells_ptr

  subroutine universe_ptr_assignment(LHS,RHS)
    class(universe_ptr), intent(out)  :: LHS
    type(universe_ptr), intent(in)    :: RHS

    if(associated(LHS % ptr)) deallocate(LHS % ptr)
    LHS % ptr => RHS % ptr
  end subroutine universe_ptr_assignment

  subroutine universe_ptr_assignment_target(LHS,RHS)
    class(universe_ptr), intent(out)        :: LHS
    class(universe), target, intent(in)     :: RHS

    if(associated(LHS % ptr)) deallocate(LHS % ptr)
    LHS % ptr => RHS
  end subroutine universe_ptr_assignment_target

  !
  ! Check whether the pointer wrapper is associated to a universe
  !
  function associated_ptr(self)result(assoc)
    class(universe_ptr), intent(in) :: self
    logical(defBool) :: assoc
    assoc = associated(self % ptr)
  end function associated_ptr

  subroutine kill(self)
    class(universe_ptr), intent(inout) :: self
    self % ptr => null()
  end subroutine kill

end module universe_class
