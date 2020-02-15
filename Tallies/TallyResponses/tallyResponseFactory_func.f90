module tallyResponseFactory_func

  use numPrecision
  use genericProcedures, only : fatalError
  use dictionary_class,  only : dictionary

  ! tallyResponse interface
  use tallyResponse_inter,    only : tallyResponse

  ! tallyResponse implementations
  use fluxResponse_class,     only : fluxResponse
  use macroResponse_class,    only : macroResponse
  use testResponse_class,     only : testResponse

  implicit none
  private

  public :: new_tallyResponse

  ! *** ADD NAME OF A NEW TALLY FILTER HERE ***!
  ! List that contains all accaptable types of tallyResponses
  ! It is printed if type was unrecognised
  ! NOTE:
  ! For now  it is necessary to adjust trailing blanks so all enteries have the same length
  character(nameLen),dimension(*),parameter :: AVALIBLE_tallyResponses = ['fluxResponse ',&
                                                                          'macroResponse']

contains

  !!
  !! Allocate new allocatable tallyResponse to a specific type
  !! If new is allocated it deallocates it
  !!
  !!
  subroutine new_tallyResponse(new,dict)
    class(tallyResponse),allocatable, intent(inout) :: new
    class(dictionary), intent(in)                   :: dict
    character(nameLen)            :: type
    character(100),parameter      :: Here = 'new_tallyResponse (tallyResponseFactory_func.f90)'

    ! Deallocate new if allocated
    if(allocated(new)) deallocate(new)

    ! Obtain string that specifies type to be built
    call dict % get(type,'type')

    ! Allocate approperiate subclass of tallyResponse
    ! *** ADD CASE STATEMENT FOR A NEW TALLY MAP BELOW ***!
    select case(type)
      case('fluxResponse')
        allocate(fluxResponse :: new)
        call new % init(dict)

      case('macroResponse')
        allocate(macroResponse :: new)
        call new % init(dict)

      case('testResponse')
        allocate(testResponse :: new)
        call new % init(dict)

     !*** NEW TALLY MAP TEMPLATE ***!
     !case('<newtallyResponseName>')
     !  allocate(<newtallyResponseName> :: new)
     !  call new % init(dict)
     !
      case default
        print *, AVALIBLE_tallyResponses
        call fatalError(Here, 'Unrecognised type of tallyResponse: ' // trim(type))

    end select

  end subroutine new_tallyResponse

end module tallyResponseFactory_func
