!!
!! Factory of tallyMap
!!
!! Builds an instance of any tallyMap from a dictionary
!!
!! Interface:
!!   new_tallyMap -> build an instance of class(tallyMap)
!!
module tallyMapFactory_func

  use numPrecision
  use universalVariables
  use genericProcedures, only : fatalError, linFind
  use dictionary_class,  only : dictionary

  ! TallyMap interface and sub-interface Factories
  use tallyMap_inter,         only : tallyMap
  use tallyMap1DFactory_func, only : new_tallyMap1D => new_tallyMap, AVALIBLE_tallyMaps1D


  ! TallyMap implementations
  use multiMap_class,         only : multiMap
  use sphericalMap_class,     only : sphericalMap
  use cylindricalMap_class,   only : cylindricalMap

  implicit none
  private

  public :: new_tallyMap


  ! *** ADD NAME OF A NEW TALLY MAP 1D HERE ***!
  ! List that contains all accaptable types of tallyMaps1D
  ! It is printed if type was unrecognised
  ! NOTE:
  ! For now  it is necessary to adjust trailing blanks so all enteries have the same length
  character(nameLen),dimension(*),parameter :: AVALIBLE_tallyMaps = [ 'multiMap      ', &
                                                                      'sphericalMap  ', &
                                                                      'cylindricalMap']


contains

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
  !!   Will return an error if type of tallyMap is not recognised
  !!
  subroutine new_tallyMap(new, dict)
    class(tallyMap),allocatable, intent(inout) :: new
    class(dictionary), intent(in)              :: dict
    character(nameLen)                         :: type
    character(100),parameter  :: Here = 'new_tallyMap (tallyMapFactory_func.f90)'

    ! Deallocate new if allocated
    if(allocated(new)) deallocate(new)

    ! Obtain string that specifies type to be built
    call dict % get(type,'type')

    ! Check if the requested type is tallyMap1D
    if (linFind(AVALIBLE_tallyMaps1D, type) /= targetNotFound) then
      call new_tallyMap1D(new, dict)

    else ! Check aginst multidimensional tallyMaps
      select case(type)
        case('multiMap')
          allocate( multiMap :: new)
          call new % init(dict)

        case('sphericalMap')
          allocate( sphericalMap :: new)
          call new % init(dict)

        case('cylindricalMap')
          allocate( cylindricalMap :: new)
          call new % init(dict)

      end select
    end if

    ! Print error if failed to build a tallyMap
    if(.not.allocated(new)) then
      print *, AVALIBLE_tallyMaps1D, AVALIBLE_tallyMaps
      call fatalError(Here,' Unrecognised tallyMap : '// trim(type))
    end if

  end subroutine new_tallyMap

end module tallyMapFactory_func
