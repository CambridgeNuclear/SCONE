module tallyMapFactory_func

  use numPrecision
  use genericProcedures, only : fatalError
  use dictionary_class,  only : dictionary

  ! TallyMap interface
  use tallyMap_inter,    only : tallyMap

  ! TallyMap implementations
  use energyMap_class,   only : energyMap
  use spaceMap_class,    only : spaceMap
  use materialMap_class, only : materialMap
  use testMap_class,     only : testMap
!  use matXsMap_class,    only : matXsMap

  implicit none
  private

  public :: new_tallyMap

  ! *** ADD NAME OF A NEW TALLY MAP HERE ***!
  ! List that contains all accaptable types of tallyMaps
  ! It is printed if type was unrecognised
  ! NOTE:
  ! For now  it is necessary to adjust trailing blanks so all enteries have the same length
  character(nameLen),dimension(*),parameter :: AVALIBLE_tallyMaps = [ 'energyMap  ',&
                                                                      'spaceMap   ',&
                                                                      'materialMap',&
                                                                      'testMap    ']

contains

  !!
  !! Allocate new allocatable tallyMap to a specific type
  !! If new is allocated it deallocates it
  !!
  !!
  subroutine new_tallyMap(new, dict)
    class(tallyMap),allocatable, intent(inout) :: new
    class(dictionary), intent(in)              :: dict
    character(nameLen)                         :: type
    character(100),parameter                   :: Here = 'new_tallyMap (tallyMapFactory_func.f90)'

    ! Deallocate new if allocated
    if(allocated(new)) deallocate(new)

    ! Obtain string that specifies type to be built
    call dict % get(type,'type')

    ! Allocate approperiate subclass of tallyMap
    ! *** ADD CASE STATEMENT FOR A NEW TALLY MAP BELOW ***!
    select case(type)
      case('energyMap')
        allocate(energyMap :: new)
        call new % init(dict)

      case('spaceMap')
        allocate(spaceMap :: new)
        call new % init(dict)

      case('materialMap')
        allocate(materialMap :: new)
        call new % init(dict)

      case('testMap')
        allocate(testMap :: new)
        call new % init(dict)

     !*** NEW TALLY MAP TEMPLATE ***!
     !case('<newTallyMapName>')
     !  allocate(<newTallyMapName> :: new)
     !  call new % init(dict)
     !
      case default
        print *, AVALIBLE_tallyMaps
        call fatalError(Here, 'Unrecognised type of tallyMap: ' // trim(type))

    end select

  end subroutine new_tallyMap

end module tallyMapFactory_func
