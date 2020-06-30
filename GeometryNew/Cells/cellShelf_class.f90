module cellShelf_class

  use numPrecision
  use genericProcedures,  only : fatalError, numToChar
  use dictionary_class,   only : dictionary
  use intMap_class,       only : intMap
  use surfaceShelf_class, only : surfaceShelf
  use cell_inter,         only : cell
  use cellFactory_func,   only : new_cell_ptr

  implicit none
  private

  !!
  !! Small, local container to store polymorphic cells in a single array
  !!
  !! Public members:
  !!   name -> Name of the cell
  !!   ptr  -> Pointer to the cell
  !!
  type :: cellBox
    character(nameLen)   :: name = ''
    class(cell), pointer :: ptr  => null()
  end type cellBox

  !!
  !! Storage space for cells defined in the geometry
  !!
  !! Sample dictionary input:
  !!   cells {
  !!     cell1 {<cell definition>}
  !!     cell2 {<cell definition>}
  !!     ...
  !!      }
  !!
  !! Private Members:
  !!   cells -> array to store pointers to polymorphic cells
  !!   idMap -> Map between cellID and corresponding index
  !!
  !! Interface:
  !!   init -> Initialise from a dictionary & surfaceShelf
  !!   cellPtr -> Get pointer to a cell given by its index
  !!   cellIdx -> Return index of a cell fivent its ID
  !!   cellID  -> Return cell ID given its index
  !!   kill    -> Return to uninitialised state
  !!
  !! NOTE: Becouse cells are stored as pointers, calling `kill` is crutial to prevent
  !!   memory leaks. TODO: Add `final` procedure here ?
  !!
  type, public :: cellShelf
    private
    type(cellBox), dimension(:), allocatable :: cells
    type(intMap)                             :: idMap
  contains
    procedure :: init
    procedure :: cellPtr
    procedure :: cellIdx
    procedure :: cellID
    procedure :: kill
  end type cellShelf

contains

  !!
  !! Load cells into the shelf
  !!
  !! Args:
  !!   dict [in] -> Dictionary with subdictionaries that contain cell definitions
  !!   surfs [inout] -> Surface shelf with user-defined surfaces
  !!
  !! Errors:
  !!   fatalError if there are clashes in surface ID
  !!
  subroutine init(self, dict, surfs)
    class(cellShelf), intent(inout)               :: self
    class(dictionary), intent(in)                 :: dict
    type(surfaceShelf), intent(inout)             :: surfs
    character(nameLen), dimension(:), allocatable :: names
    integer(shortInt)                             :: i, id, idx
    integer(shortInt), parameter :: NOT_PRESENT = -7
    character(100), parameter :: Here = 'init (cellShelf_class.f90)'

    ! Get all keys for subdictionaries
    call dict % keys(names, 'dict')

    ! Allocate space
    allocate (self % cells(size(names)))

    ! Build cells
    do i = 1, size(names)
      self % cells(i) % name = names(i)
      self % cells(i) % ptr => new_cell_ptr(dict % getDictPtr(names(i)), surfs)
      id = self % cells(i) % ptr % id()

      ! Add ID to the map detecting any conflicts
      idx = self % idMap % getOrDefault(id, NOT_PRESENT)
      if (idx /= NOT_PRESENT) then
        call fatalError(Here,'Cells '//trim(names(i))// ' & '//&
                              trim(self % cells(idx) % name)//&
                             ' have the same ID: '//numToChar(id))

      else
        call self % idMap % add(id, i)

      end if
    end do

  end subroutine init

  !!
  !! Return pointer to the cell indicated by index
  !!
  !! Args:
  !!   idx [in] -> Index of the cell
  !!
  !! Result:
  !!   Pointer to the cell under index idx
  !!
  !! Error:
  !!   fatalError is idx does not correspond to a cell (is out-of-bounds)
  !!
  function cellPtr(self, idx) result (ptr)
    class(cellShelf), intent(in)  :: self
    integer(shortInt), intent(in) :: idx
    class(cell), pointer          :: ptr
    character(100), parameter :: Here = 'cellPtr (cellShelf_class.f90)'

    ! Catch invalid idx
    if (idx < 1 .or. idx > size(self % cells)) then
       call fatalError(Here, 'Requested index: '//numToChar(idx)//' is not valid. Must be between &
                             &1 and '//numToChar(size(self % cells)))
    end if

    ! Return pointer
    ptr => self % cells(idx) % ptr

  end function cellPtr

  !!
  !! Return IDX of a cell with ID
  !!
  !! Args:
  !!   id [in] -> Id of the cell
  !!
  !! Result:
  !!   Index of the cell with ID
  !!
  !! Error:
  !!   fatalError if there is no cell with ID
  !!
  function cellIdx(self, id) result(idx)
    class(cellShelf), intent(in)    :: self
    integer(shortInt), intent(in)   :: id
    integer(shortInt)               :: idx
    integer(shortInt), parameter :: NOT_PRESENT = -7
    character(100), parameter :: Here = 'cellIdx (cellShelf_class.f90)'

    idx = self % idMap % getOrDefault(id, NOT_PRESENT)

    if (idx == NOT_PRESENT) then
      call fatalError(Here, 'There is no cell with ID: '//numToChar(id))
    end if

  end function cellIdx

  !!
  !! Return ID of the cell with index
  !!
  !! Args:
  !!   idx [in] -> Index of the cell
  !!
  !! Result:
  !!   ID of the cell under index
  !!
  !! Error:
  !!   fatalError is idx does not correspond to a cell (is out-of-bounds)
  !!
  function cellId(self, idx) result(id)
    class(cellShelf), intent(in)  :: self
    integer(shortInt), intent(in) :: idx
    integer(shortInt)             :: id
    character(100), parameter :: Here = 'cellID (cellShelf_class.f90)'

    ! Catch invalid idx
    if (idx < 1 .or. idx > size(self % cells)) then
       call fatalError(Here, 'Requested index: '//numToChar(idx)//' is not valid. Must be between &
                             &1 and '//numToChar(size(self % cells)))
    end if

    id = self % cells(idx) % ptr % id()

  end function cellId

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(cellShelf), intent(inout) :: self
    integer(shortInt)                  :: i

    if (allocated(self % cells)) then
      do i = 1, size(self % cells)
        call self % cells(i) % ptr % kill()
      end do

      deallocate(self % cells)
    end if

    call self % idMap % kill()

  end subroutine kill

end module cellShelf_class
