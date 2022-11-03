module cellShelf_class

  use numPrecision
  use universalVariables, only : OUTSIDE_MAT
  use genericProcedures,  only : fatalError, numToChar
  use dictionary_class,   only : dictionary
  use intMap_class,       only : intMap
  use charMap_class,      only : charMap
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
  !! Also stores information about the content of the cells in terms of uniID or matIdx.
  !! uniIDs are stored as negative numbers and matIdx as positive.
  !!
  !! Sample dictionary input:
  !!   cells {
  !!     cell1 {<cell definition>; filltype uni; universe 17;}
  !!     cell2 {<cell definition>; filltype mat; material void;}
  !!     cell3 {<cell definition>; filltype outside; }
  !!     ...
  !!      }
  !!
  !! Private Members:
  !!   cells -> array to store pointers to polymorphic cells
  !!   idMap -> Map between cellID and corresponding index
  !!   fillMap -> Map between cell IDX (index) and its fillng (-uniID or matIdx)
  !!
  !! Interface:
  !!   init    -> Initialise from a dictionary & surfaceShelf
  !!   getPtr  -> Get pointer to a cell given by its index
  !!   getIdx  -> Return index of a cell given its ID
  !!   getID   -> Return cell ID given its index
  !!   getFill -> Return content of the cell. If -ve it is universe ID. If +ve it is matIdx.
  !!   getSize -> Return the number of cells (max cellIdx)
  !!   kill    -> Return to uninitialised state
  !!
  !! NOTE: Because cells are stored as pointers, calling `kill` is crucial to prevent
  !!   memory leaks. TODO: Add `final` procedure here ?
  !!
  type, public :: cellShelf
    private
    type(cellBox), dimension(:), allocatable :: cells
    type(intMap)                             :: idMap
    type(intMap)                             :: fillMap
  contains
    procedure :: init
    procedure :: getPtr
    procedure :: getIdx
    procedure :: getID
    procedure :: getFill
    procedure :: getSize
    procedure :: kill
  end type cellShelf

contains

  !!
  !! Load cells into the shelf
  !!
  !! Args:
  !!   dict [in] -> Dictionary with subdictionaries that contain cell definitions
  !!   surfs [inout] -> Surface shelf with user-defined surfaces
  !!   mats [in] -> Map of material names to matIdx
  !!
  !! Errors:
  !!   fatalError if there are clashes in cell ID
  !!
  subroutine init(self, dict, surfs, mats)
    class(cellShelf), intent(inout)               :: self
    class(dictionary), intent(in)                 :: dict
    type(surfaceShelf), intent(inout)             :: surfs
    type(charMap), intent(in)                     :: mats
    character(nameLen), dimension(:), allocatable :: names
    class(dictionary), pointer                    :: tempDict
    character(nameLen)                            :: filling, matName
    integer(shortInt)                             :: i, id, idx, fill
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

      ! Load cell content
      tempDict => dict % getDictPtr(names(i))
      call tempDict % get(filling, 'filltype')

      select case (filling)
        case ('outside')
          call self % fillMap % add(i, OUTSIDE_MAT)

        case ('mat')
          call tempDict % get(matName, 'material')
          fill = mats % getOrDefault(matName, NOT_PRESENT)

          if (fill == NOT_PRESENT) then
            call fatalError(Here, 'Material with name '//trim(matName)//' was not found.')
          end if
          call self % fillMap % add(i, fill)

        case ('uni')
          call tempDict % get(fill, 'universe')

          if (fill <= 0) then
            call fatalError(Here, 'Universe ID must be +ve. Is :'//numToChar(fill))
          end if
          call self % fillMap % add(i, -fill)

        case default
          call fatalError(Here, 'Unknown type of cell filling: '//trim(filling))

      end select
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
  function getPtr(self, idx) result (ptr)
    class(cellShelf), intent(in)  :: self
    integer(shortInt), intent(in) :: idx
    class(cell), pointer          :: ptr
    character(100), parameter :: Here = 'getPtr (cellShelf_class.f90)'

    ! Catch invalid idx
    if (idx < 1 .or. idx > size(self % cells)) then
       call fatalError(Here, 'Requested index: '//numToChar(idx)//' is not valid. Must be between &
                             &1 and '//numToChar(size(self % cells)))
    end if

    ! Return pointer
    ptr => self % cells(idx) % ptr

  end function getPtr

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
  function getIdx(self, id) result(idx)
    class(cellShelf), intent(in)    :: self
    integer(shortInt), intent(in)   :: id
    integer(shortInt)               :: idx
    integer(shortInt), parameter :: NOT_PRESENT = -7
    character(100), parameter :: Here = 'getIdx (cellShelf_class.f90)'

    idx = self % idMap % getOrDefault(id, NOT_PRESENT)

    if (idx == NOT_PRESENT) then
      call fatalError(Here, 'There is no cell with ID: '//numToChar(id))
    end if

  end function getIdx

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
  function getId(self, idx) result(id)
    class(cellShelf), intent(in)  :: self
    integer(shortInt), intent(in) :: idx
    integer(shortInt)             :: id
    character(100), parameter :: Here = 'getId (cellShelf_class.f90)'

    ! Catch invalid idx
    if (idx < 1 .or. idx > size(self % cells)) then
       call fatalError(Here, 'Requested index: '//numToChar(idx)//' is not valid. Must be between &
                             &1 and '//numToChar(size(self % cells)))
    end if

    id = self % cells(idx) % ptr % id()

  end function getId

  !!
  !! Obtain contents of the cell given by index
  !!
  !! Args:
  !!   idx [in] -> Index of the cell
  !!
  !! Result:
  !!   Content of the cell. Negative (-ve) values represent universe ID (-uniId).
  !!   Positive values are material indexes (matIdx).
  !!
  function getFill(self, idx) result(fill)
    class(cellShelf), intent(in)    :: self
    integer(shortInt), intent(in)   :: idx
    integer(shortInt)               :: fill
    integer(shortInt), parameter :: NOT_PRESENT = -7
    character(100), parameter :: Here = 'getFill (cellShelf_class.f90)'

    ! Catch invalid idx
    if (idx < 1 .or. idx > size(self % cells)) then
       call fatalError(Here, 'Requested index: '//numToChar(idx)//' is not valid. Must be between &
                             &1 and '//numToChar(size(self % cells)))
    end if

    ! Obtain fill
    fill = self % fillMap % get(idx)

  end function getFill

  !!
  !! Return size of the shelf
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   Number of cells on the shelf
  !!
  elemental function getSize(self) result(N)
    class(cellShelf), intent(in) :: self
    integer(shortInt)            :: N

    N = size(self % cells)

  end function getSize

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
    call self % fillMap % kill()

  end subroutine kill

end module cellShelf_class
