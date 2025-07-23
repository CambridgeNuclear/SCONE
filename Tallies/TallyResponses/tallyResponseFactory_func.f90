module tallyResponseFactory_func

  use numPrecision
  use genericProcedures, only : fatalError
  use dictionary_class,  only : dictionary

  ! tallyResponse interface
  use tallyResponse_inter,    only : tallyResponse

  ! tallyResponse implementations
  use fluxResponse_class,     only : fluxResponse
  use macroResponse_class,    only : macroResponse
  use microResponse_class,    only : microResponse
  use weightResponse_class,   only : weightResponse
  use invSpeedResponse_class, only : invSpeedResponse
  use testResponse_class,     only : testResponse

  implicit none
  private

  public :: new_tallyResponse

  ! List that contains all accaptable types of tallyResponses
  ! It is printed if type was unrecognised
  ! NOTE:
  ! For now  it is necessary to adjust trailing blanks so all enteries have the same length
  character(nameLen),dimension(*),parameter :: AVALIBLE_tallyResponses = ['fluxResponse    ',&
                                                                          'macroResponse   ',&
                                                                          'microResponse   ',&
                                                                          'weightResponse  ',&
                                                                          'invSpeedResponse']

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
    select case(type)
      case('fluxResponse')
        allocate(fluxResponse :: new)

      case('macroResponse')
        allocate(macroResponse :: new)

      case('microResponse')
        allocate(microResponse :: new)

      case('weightResponse')
        allocate(weightResponse :: new)

      case('invSpeedResponse')
        allocate(invSpeedResponse :: new)

      case('testResponse')
        allocate(testResponse :: new)

      case default
        print *, AVALIBLE_tallyResponses
        call fatalError(Here, 'Unrecognised type of tallyResponse: ' // trim(type))

    end select

  ! Initialise new response
  call new % init(dict)

  end subroutine new_tallyResponse

end module tallyResponseFactory_func
