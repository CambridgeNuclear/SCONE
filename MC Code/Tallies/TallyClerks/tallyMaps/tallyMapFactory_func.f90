module tallyMapFactory_func

  use numPrecision
  use genericProcedures, only : fatalError
  use dictionary_class,  only : dictionary

  ! TallyMap interface
  use tallyMap_inter,    only : tallyMap

  ! TallyMap implementations
  use energyMap_class,   only : energyMap
  use spaceMap_class,    only : spaceMap

  implicit none
  private

  public :: new_tallyMap
  public :: new_tallyMap_ptr

  ! *** ADD NAME OF A NEW TALLY MAP HERE ***!
  ! List that contains all accaptable types of tallyMaps
  ! It is printed if type was unrecognised
  ! NOTE:
  ! For now  it is necessary to adjust trailing blanks so all enteries have the same length
  character(nameLen),dimension(*),parameter :: AVALIBLE_tallyMaps = [ 'energyMap',&
                                                                      'spaceMap ' ]

contains

  !!
  !! Returns allocatable tallyMap form dictionary
  !!
  function new_tallyMap(dict) result(new)
    class(dictionary), intent(in) :: dict
    class(tallyMap),allocatable   :: new
    character(nameLen)            :: type
    character(100),parameter      :: Here = 'new_tallyMap (tallyMapFactory_func.f90)'

    ! Obtain string that specifies type to be built
    call dict % get(type,'type')

    ! Allocate approperiate subclass of tallyMap
    ! *** ADD CASE STATEMENT FOR A NEW TALLY MAP BELOW ***!
    select case(type)
      case('energyMap')
        allocate(new, source = energyMap(dict) )

      case('spaceMap')
        allocate(new, source = spaceMap(dict))

     !*** NEW TALLY MAP TEMPLATE ***!
     !case('<newTallyMapName>')
     !  allocate(new, source = <newTallyMapName>(dict) )
     !
      case default
        print *, AVALIBLE_tallyMaps
        call fatalError(Here, 'Unrecognised type of tallyMap: ' // trim(type))

    end select


  end function new_tallyMap

  !!
  !! Returns pointer to allocated tallyMap from dictionary
  !!
  function new_tallyMap_ptr(dict) result(new)
    class(dictionary), intent(in) :: dict
    class(tallyMap),pointer        :: new

    ! Allocate pointer and copy data from local allocatable
    allocate( new, source = new_tallyMap(dict) )

  end function new_tallyMap_ptr
    
end module tallyMapFactory_func
