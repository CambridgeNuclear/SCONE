module universeFactory_func

  use numPrecision
  use genericProcedures,  only : fatalError
  use dictionary_class,   only : dictionary
  use maps_class,         only : intMap

  ! Surface and Cell Shelf
  use surface_inter,      only : surfaceShelf
  use cell_class,         only : cellShelf

  ! Universe interface
  use universe_inter,     only : universe

  ! Universe implementations
  use cellUniverse_class, only : cellUniverse
  use pinUniverse_class,  only : pinUniverse
  use latUniverse_class,  only : latUniverse

  !*** STAYS HERE ONLY PROVISIONALLY
  use nuclearData_inter,  only : nuclearData

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
  character(nameLen),dimension(*),parameter :: AVALIBLE_universes = [ 'cellUniverse',&
                                                                      'pinUniverse ',&
                                                                      'latUniverse ']

contains

  !!
  !! Returns allocatable universe form dictionary and surface & cell Shelfs
  !! Return fillVector as well. +ve entries are material Idxs, -ve are fill universes Ids
  !!
  function new_universe(fillVector, dict, cShelf, sShelf, cellFillMap ,materials ) result(new)
    integer(shortInt),dimension(:),allocatable,intent(out) :: fillVector
    class(dictionary), intent(in)                          :: dict
    type(cellShelf), intent(inout)                         :: cShelf
    type(surfaceShelf), intent(inout)                      :: sShelf
    type(intMap), intent(in)                               :: cellFillMap
    class(nuclearData), intent(in)                         :: materials
    class(universe),allocatable                            :: new
    character(nameLen)                                     :: type
    character(100),parameter          :: Here = 'new_universe (universeFactory_func.f90)'

    ! Obtain string that specifies type to be built
    call dict % get(type,'type')

    ! Allocate approperiate subclass of universe
    ! *** ADD CASE STATEMENT FOR A NEW SURFACE BELOW ***!
    select case(type)
      case('cellUniverse')
        allocate(new, source = cellUniverse(fillVector, dict, cShelf, sShelf, cellFillMap,materials))

      case('pinUniverse')
        allocate(new, source = pinUniverse(fillVector, dict, cShelf, sShelf, cellFillMap,materials))

      case('latUniverse')
        allocate(new, source = latUniverse(fillVector, dict, cShelf, sShelf, cellFillMap,materials))

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
  !! Return fillVector as well. +ve entries are material Idxs, -ve are fill universes Ids
  !!
  function new_universe_ptr(fillVector,dict, cShelf, sShelf, cellFillMap ,materials) result(new)
    integer(shortInt),dimension(:),allocatable,intent(out) :: fillVector
    class(dictionary), intent(in)                          :: dict
    type(cellShelf), intent(inout)                         :: cShelf
    type(surfaceShelf), intent(inout)                      :: sShelf
    type(intMap), intent(in)                               :: cellFillMap
    class(nuclearData), intent(in)                         :: materials
    class(universe),pointer                                :: new

    ! Allocate pointer and copy data from local allocatable
    allocate( new, source = new_universe(fillVector ,dict, cShelf, sShelf, cellFillMap, materials) )

  end function new_universe_ptr


end module universeFactory_func
