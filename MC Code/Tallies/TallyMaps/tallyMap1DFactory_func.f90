!!
!! Factory of tallyMap1D
!!
!! Build an instance of a tallyMap1D from a dictionary. Can build it as class(tallyMap1D) or
!! as class(tallyMap)
!!
!! Public Members:
!!   AVALIBLE_tallyMaps1D -> PARAMETER nameLen array of names of avalible talyMaps1D
!!
!! Interface:
!!   new_tallyMap1D -> builds an instance as a class(tallyMap1D)
!!   new_tallyMap   -> builds an instance as a class(tallyMap)
!!
!! NOTE:
!!   This factory is used by general tallyMap factory
!!
module tallyMap1DFactory_func

  use numPrecision
  use genericProcedures, only : fatalError
  use dictionary_class,  only : dictionary

  ! TallyMap interface
  use tallyMap_inter,      only : tallyMap
  use tallyMap1D_inter,    only : tallyMap1D

  ! TallyMap implementations
  use energyMap_class,   only : energyMap
  use spaceMap_class,    only : spaceMap
  use materialMap_class, only : materialMap
  use testMap_class,     only : testMap
!  use matXsMap_class,    only : matXsMap

  implicit none
  private

  public :: new_tallyMap1D
  public :: new_tallyMap


  ! *** ADD NAME OF A NEW TALLY MAP 1D HERE ***!
  ! List that contains all accaptable types of tallyMaps1D
  ! It is printed if type was unrecognised
  ! NOTE:
  ! For now  it is necessary to adjust trailing blanks so all enteries have the same length
  character(nameLen),dimension(*),parameter, public :: AVALIBLE_tallyMaps1D = [ 'energyMap  ',&
                                                                                'spaceMap   ',&
                                                                                'materialMap',&
                                                                                'testMap    ']

contains

  !!
  !! Allocates new allocatable tallyMap to a specific type as class(tallyMap1D)
  !! If new is allocated it deallocates it
  !!
  !! Args:
  !!   new [inout] -> an allocatable class(tallyMap1D), will be allocated on exit. Any
  !!                  existing content will be deallocated (NO kill subroutine called !!)
  !!   dict [in]   -> dictionary with the settings
  !!
  !! Errors:
  !!   Will return an error if type of tallyMap1D is not recognised
  !!
  subroutine new_tallyMap1D(new, dict)
    class(tallyMap1D),allocatable, intent(inout) :: new
    class(dictionary), intent(in)                :: dict
    character(nameLen)                           :: type
    character(100),parameter  :: Here = 'new_tallyMap1D (tallyMap1DFactory_func.f90)'

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
    end select

    ! Print error if failed to allocate
    if(.not.allocated(new)) then
      print *, AVALIBLE_tallyMaps1D
      call fatalError(Here,'Unrecognised type of tallyMap1D : ' // trim(type))
    end if

  end subroutine new_tallyMap1D

  !!
  !! Allocates new allocatable tallyMap to a specific type as class(tallyMap)
  !! If new is allocated it deallocates it
  !!
  !! Args:
  !!   new [inout] -> an allocatable class(tallyMap), will be allocated on exit. Any
  !!                  existing content will be deallocated (NO kill subroutine called !!)
  !!   dict [in]   -> dictionary with the settings
  !!
  !! Errors:
  !!   Any error given by new_tallyMap1D
  !!
  subroutine new_tallyMap(new, dict)
    class(tallyMap), allocatable, intent(inout) :: new
    class(dictionary), intent(in)               :: dict
    class(tallyMap1D), allocatable              :: temp

    ! Deallocate if allocated
    if(allocated(new)) deallocate(new)

    ! Build an instance
    call new_tallyMap1D(temp, dict)
    call move_alloc(temp, new)

  end subroutine new_tallyMap

end module tallyMap1DFactory_func
