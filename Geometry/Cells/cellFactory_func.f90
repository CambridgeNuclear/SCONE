module cellFactory_func

  use numPrecision
  use genericProcedures,  only : fatalError
  use dictionary_class,   only : dictionary
  use surfaceShelf_class, only : surfaceShelf

  ! Cell interface
  use cell_inter,         only : cell

  ! Cells
  use simpleCell_class,   only : simpleCell
  use unionCell_class,    only : unionCell

  implicit none
  private

  ! List that contains acceptable types of cells
  ! NOTE: It is necessary to adjust trailing blanks so all entries have the same length
  character(nameLen), dimension(*), parameter :: AVAILABLE_CELL = ['simpleCell', &
                                                                   'unionCell ']

  ! Public Interface
  public :: new_cell_ptr
  public :: new_cell

contains

  !!
  !! Return a pointer to a new instance of an allocated cell
  !!
  !! Args:
  !!   dict [in] -> Dictionary with cell definition
  !!   surfs [inout] -> Surface shelf with surfaces already defined
  !!
  !! Result:
  !!   class(cell) pointer to an allocated instance of a cell
  !!
  !! Errors:
  !!   fatalError if type of cell is unknown
  !!
  function new_cell_ptr(dict, surfs) result(new)
    class(dictionary), intent(in)     :: dict
    type(surfaceShelf), intent(inout) :: surfs
    class(cell), pointer              :: new
    character(nameLen)                :: type
    character(100), parameter :: Here = 'new_cell_ptr (cellFactory_func.f90)'

    ! Obtain type of the cell
    call dict % get(type, 'type')

    ! Allocate approperiate cell
    select case (type)
      case ('simpleCell')
        allocate(simpleCell :: new)
      
      case ('unionCell')
        allocate(unionCell :: new)

      case default
        print '(A)', 'AVAILABLE CELLS: '
        print '(A)', AVAILABLE_CELL
        call fatalError(Here, 'Unrecognised type of a cell: '//trim(type))

    end select

    ! Initialise cell
    call new % init(dict, surfs)

  end function new_cell_ptr

  !!
  !! Allocate an allocatable cell
  !!
  !! Args:
  !!   new [out]     -> Cell to be allocated
  !!   dict [in]     -> Dictionary with the cell definition
  !!   surfs [inout] -> Surface shelf with user-defined surfaces
  !!
  !! Errors:
  !!   fatalError if type of cell is unrecognised
  !!
  subroutine new_cell(new, dict, surfs)
    class(cell), allocatable, intent(out) :: new
    class(dictionary), intent(in)         :: dict
    type(surfaceShelf), intent(inout)     :: surfs
    class(cell), pointer                  :: temp

    temp => new_cell_ptr(dict, surfs)
    allocate (new, source = temp)
    call temp % kill()
    deallocate (temp)

  end subroutine new_cell


end module cellFactory_func
