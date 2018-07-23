module universeFactory_func

  use numPrecision
  use genericProcedures,  only : fatalError
  use dictionary_class,   only : dictionary

  ! Surface and Cell Shelf
  use surface_inter,      only : surfaceShelf
  use cell_class,         only : cellShelf

  ! Universe interface
  use universe_inter,     only : universe

  ! Universe implementations
  use cellUniverse_class, only : cellUniverse

  implicit none
  private

  ! Public functions to build new universe
  public :: new_universe
  public :: new_universe_ptr


  ! *** ADD NAME OF A NEW UNIVERSE HERE ***!
  ! List that contains all accaptable types of universes
  ! It is printed if type was unrecognised
  ! NOTE:
  ! For now  it is necessary to adjust trailing blanks so all enteries have the same length
  character(nameLen),dimension(*),parameter :: AVALIBLE_universes = [ 'cellUniverse']

contains

  !!
  !! Returns allocatable universe form dictionary and surface & cell Shelfs
  !! Return fillVector as well. +ve entries are fill cells, -ve are fill universes
  !!
  function new_universe(fillVector, dict, cShelf, sShelf ) result(new)
    integer(shortInt),dimension(:),allocatable,intent(out) :: fillVector
    class(dictionary), intent(in)                          :: dict
    type(cellShelf), intent(inout)                         :: cShelf
    type(surfaceShelf), intent(inout)                      :: sShelf
    class(universe),allocatable                            :: new
    character(nameLen)                                     :: type
    character(100),parameter          :: Here = 'new_universe (universeFactory_func.f90)'

    ! Obtain string that specifies type to be built
    call dict % get(type,'type')

    ! Allocate approperiate subclass of universe
    ! *** ADD CASE STATEMENT FOR A NEW SURFACE BELOW ***!
    select case(type)
      case('cellUniverse')
        allocate(new, source = cellUniverse(fillVector, dict, cShelf, sShelf)

     !*** NEW SURFACE TEMPLATE ***!
     !case('<newSUrfaceName>')
     !  allocate(new, source = <newSurfaceName>(dict) )
     !
      case default
        print *, AVALIBLE_universes
        call fatalError(Here, 'Unrecognised type of universe: ' // trim(type))

    end select

  end function new_universe

  !!
  !! Returns pointer to allocated universe from dictionary and surface & cell Shelfs
  !! Return fillVector as well. +ve entries are fill cells, -ve are fill universes
  !!
  function new_universe_ptr(fillVector, cShelf, sShelf) result(new)
    integer(shortInt),dimension(:),allocatable,intent(out) :: fillVector
    class(dictionary), intent(in)     :: dict
    type(cellShelf), intent(inout)    :: cShelf
    type(surfaceShelf), intent(inout) :: sShelf
    class(universe),pointer           :: new

    ! Allocate pointer and copy data from local allocatable
    allocate( new, source = new_surface(fillVector ,dict, cShelf, sShelf) )

  end function new_universe_ptr


end module universeFactory_func
