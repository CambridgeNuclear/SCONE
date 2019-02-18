!!
!! This module works as an object(class) with a single instance [Singelton]
!!  It purpose is to store diffrent definition of the energy grids including
!!  a number of predefined energy grids.
!!
!! Sample energyGrid definition dictionaries:
!!
!!  structGridName {
!!    type lin|log;
!!    min  <real>;
!!    max  <real>;
!!    N    <int>;
!!  }
!!
!!  unstructGridName {
!!    type unstruct;
!!    bins (<bound1>, <bound2>, ... );
!!  }
!!
!!  Following functions exist in the interface:
!!    define_energyGrid(name, dict)  -> allows to create new global named energy grid
!!    get_energyGrid(grid, name)     -> returns a "grid" objects with definition indicated by name
!!    kill_energyGrid()              -> cleans all defined energyGrids from memory
!!
!!  To add a new pre-defined energy grid it is necessary to:
!!    1) add new pre-defined energy name to PRE_DEF_NAMES array
!!    2) write new parameter array with sorted (increasing) bin boundaries in preDefEnergyGrids.f90
!!    3) add a new case to select case in get_energyGrid
!!
module energyGridRegistry_mod

  use numPrecision
  use genericProcedures, only : fatalError, charCmp, isSorted
  use dictionary_class,  only : dictionary
  use grid_class,        only : grid
  use charMap_class,     only : charMap

  implicit none
  private

  !! Map from grid name to grid index
  type(charMap) :: nameMap

  !! Stored prototypes of named energy grids and their names
  type(grid),dimension(:),allocatable :: eGrids
  integer(shortInt)                   :: eGrid_top = 1

  !! Names of predefined energy grids
  character(*),dimension(*),parameter :: PRE_DEF_NAMES = ['wims69 ',&
                                                          'wims172']

  !! Misc. parameters
  real(defReal),parameter     :: GROWTH_RATIO = 1.2_defReal
  integer(shortInt),parameter :: MIN_SIZE = 5

contains

  !!
  !! Retrieve a named energyGrid
  !!  If logical err is provided fatalError is suppresed
  !!  err is set to TRUE if subroutine FAILS to obtain energy grid
  !!
  subroutine get_energyGrid(eGrid, name, err)
    type(grid), intent(inout)             :: eGrid
    character(nameLen), intent(in)        :: name
    logical(defBool),intent(out),optional :: err
    integer(shortInt)                     :: idx
    character(100), parameter :: Here = 'get_energyGrid (energyGridRegistry_mod.f90)'

    ! Try to find the grid in the user-defined grids
    idx = nameMap % getOrDefault(name,-17)

    if( idx /= -17) then
      eGrid = eGrids(idx)
      return

    end if

    ! Try to find grid in the pre-defined structures
    ! ADD A NEW PRE_DEFINED STRUCTURE HERE
    select case(name)

    end select

    if (present(err)) then
      err = .true.

    else
      call fatalError(Here,'Grid '//name//' is undefined!')

    end if

  end subroutine get_energyGrid

  !!
  !! Defines new named energy grid
  !! Requires new grid name and dictionary with its definition
  !!
  subroutine define_energyGrid(name, dict)
    character(nameLen), intent(in)      :: name
    class(dictionary), intent(in)       :: dict
    integer(shortInt)                   :: idx, N, N_new
    type(grid),dimension(:),allocatable :: tempGrid
    character(100), parameter :: Here = 'define_energyGrid (energyGridRegistry_mod.f90)'

    ! Verify that new name does not clash with a predefined structure
    if(any(charCmp(PRE_DEF_NAMES, name))) then
      call fatalError(Here, name //' clashes with a predefined energy structure')
    end if

    ! Verify that new name was not yet defined
    if(-17 /= nameMap % getOrDefault(name, -17)) then
      call fatalError(Here, name //' energy structure was already defined')
    end if

    ! Set size of the array for storing grids
    if( allocated(eGrids)) then
      N = size(eGrids)
    else
      N = 0
    end if

    ! Grow array if needed
    if( eGrid_top >= N ) then
      ! Calculate new size
      N_new = max(MIN_SIZE, int(N * GROWTH_RATIO, shortInt))

      allocate(tempGrid(N_new))
      tempGrid(1:N) = eGrids(1:N)
      call move_alloc(tempGrid, eGrids)

    end if

    ! Build new instance of an energy grid
    call new_energyGrid(eGrids(eGrid_top), name, dict)
    call nameMap % add(name, eGrid_top)

    eGrid_top = eGrid_top + 1

  end subroutine define_energyGrid

  !!
  !! Private factory for energy structures from dictionary
  !!
  subroutine new_energyGrid(eGrid, name, dict)
    type(grid), intent(inout)               :: eGrid
    character(nameLen), intent(in)          :: name
    class(dictionary), intent(in)           :: dict
    character(nameLen)                      :: type
    integer(shortInt)                       :: N
    real(defReal)                           :: mini, maxi
    real(defReal),dimension(:), allocatable :: bins
    character(100),parameter :: Here ='new_energyGrid (energyGridRegistry_mod.f90)'

    ! Read grid type
    call dict % get(type,'grid')

    select case(type)
      case('lin', 'log')
        ! Get data
        call dict % get(N,'size')
        call dict % get(mini,'min')
        call dict % get(maxi,'max')

        ! Verify data
        if(any( [mini, maxi] < ZERO)) call fatalError(Here,' Energy grid '//name//' contains -ve energies')
        if( N <= 0 ) call fatalError(Here,'Non-positive size of the grid '//name)

        ! Build brid
        call eGrid % init(mini, maxi, N, type)

      case('unstruct')
        ! Get data
        call dict % get(bins,'bins')

        ! Verify data
        if(any(bins < ZERO)) call fatalError(Here,' Energy grid '//name//' contains -ve energies')
        if(.not.isSorted(bins)) call fatalError(Here,'Bins boundaries for grid '//name//' are not sorted')

        ! Build grid
        call eGrid % init(bins)

      case default
        call fatalError(Here,'Unknown energy grid structure: '//type//' must be in {lin;log;unstruct}')

    end select
  end subroutine new_energyGrid

    
end module energyGridRegistry_mod
