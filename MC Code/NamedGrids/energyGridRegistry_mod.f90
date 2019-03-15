!!
!! This module works as an object(class) with a single instance [Singelton]
!!  It purpose is to store diffrent definition of the energy grids including
!!  a number of predefined energy grids.
!!
!! Sample energyGrid definition dictionaries:
!!
!!  structGridName {
!!    gird lin|log;
!!    min  <real>;
!!    max  <real>;
!!    size <int>;
!!  }
!!
!!  unstructGridName {
!!    grid unstruct;
!!    bins (<bound1>, <bound2>, ... );
!!  }
!!
!!  Following functions exist in the interface:
!!    define_energyGrid(name, dict)   -> allows to create new global named energy grid
!!    define_multipleEnergyGrids(dict)-> allows to define many energyGrids defined as subdictionaries
!!                                       of the provided dictionary
!!    get_energyGrid(grid, name)      -> returns a "grid" objects with definition indicated by name
!!    kill_energyGrids()              -> cleans all defined energyGrids from memory
!!
!!  To add a new pre-defined energy grid it is necessary to:
!!    1) add new pre-defined energy name to PRE_DEF_NAMES array
!!    2) write new parameter array with sorted (increasing) bin boundaries in preDefEnergyGrids.f90
!!    3) add a new case to select case in get_energyGrid
!!
module energyGridRegistry_mod

  use numPrecision
  use preDefEnergyGrids
  use genericProcedures, only : fatalError, charCmp, isDescending
  use dictionary_class,  only : dictionary
  use energyGrid_class,  only : energyGrid
  use charMap_class,     only : charMap

  implicit none
  private

  !! Module public interface
  public :: define_energyGrid
  public :: define_multipleEnergyGrids
  public :: get_energyGrid
  public :: kill_energyGrids


  !! Map from energyGrid name to grid index
  type(charMap) :: nameMap

  !! Stored prototypes of named energy energyGrids and their names
  type(energyGrid),dimension(:),allocatable :: eGrids
  integer(shortInt)                         :: eGrid_top = 1

  !! Names of predefined energy grids
  character(*),dimension(*),parameter :: PRE_DEF_NAMES = ['wims69 ',&
                                                          'wims172',&
                                                          'casmo40']

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
    type(energyGrid), intent(inout)       :: eGrid
    character(nameLen), intent(in)        :: name
    logical(defBool),intent(out),optional :: err
    integer(shortInt)                     :: idx
    character(100), parameter :: Here = 'get_energyGrid (energyGridRegistry_mod.f90)'

    ! Try to find the grid in the user-defined grids
    idx = nameMap % getOrDefault(name,-17)
    if(present(err)) err = .false.

    if( idx /= -17) then
      eGrid = eGrids(idx)
      return

    end if

    ! Try to find grid in the pre-defined structures
    ! ADD A NEW PRE_DEFINED STRUCTURE HERE
    select case(name)
      case('wims69')
        call eGrid % init(wims69)

      case('wims172')
        call eGrid % init(wims172)

      case('casmo40')
        call eGrid % init(casmo40)

      case default
        if (present(err)) then
          err = .true.
        else
          call fatalError(Here,'Grid '//name//' is undefined!')
        end if
    end select



  end subroutine get_energyGrid

  !!
  !! Defines new named energy grid
  !! Requires new grid name and dictionary with its definition
  !!
  subroutine define_energyGrid(name, dict)
    character(nameLen), intent(in)            :: name
    class(dictionary), intent(in)             :: dict
    integer(shortInt)                         :: N, N_new
    type(energyGrid),dimension(:),allocatable :: tempGrid
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
  !! Assumes that all subdictionaries of the dictionary contain energy grids and
  !! defines them all.
  !!
  subroutine define_multipleEnergyGrids(dict)
    class(dictionary), intent(in)                :: dict
    character(nameLen), dimension(:),allocatable :: names
    integer(shortInt)                            :: i

    ! Load all dictionary names
    call dict % keys(names,'dict')

    ! Load all dictionaries
    do i=1,size(names)
      call define_energyGrid(names(i), dict % getDictPtr(names(i)))

    end do

  end subroutine define_multipleEnergyGrids

  !!
  !! Kill all definitions of energy grids
  !! Release the memory. Return to initial state
  !!
  subroutine kill_energyGrids()

    ! Grids name map
    call nameMap % kill()

    ! Grids prototypes
    if(allocated(eGrids)) then
      call eGrids % kill()
      deallocate(eGrids)
      eGrid_top = 1
    end if

  end subroutine kill_energyGrids


  !!
  !! Private factory for energy structures from dictionary
  !!
  subroutine new_energyGrid(eGrid, name, dict)
    type(energyGrid), intent(inout)               :: eGrid
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
        if(.not.isDescending(bins)) call fatalError(Here,'Bins boundaries for grid '//name//' are not descending')

        ! Build grid
        call eGrid % init(bins)

      case default
        call fatalError(Here,'Unknown energy grid structure: '//type//' must be in {lin;log;unstruct}')

    end select
  end subroutine new_energyGrid

    
end module energyGridRegistry_mod
