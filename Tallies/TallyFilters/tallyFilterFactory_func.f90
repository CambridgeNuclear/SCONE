module tallyFilterFactory_func

  use numPrecision
  use genericProcedures, only : fatalError
  use dictionary_class,  only : dictionary

  ! tallyFilter interface
  use tallyFilter_inter,    only : tallyFilter

  ! tallyFilter implementations
  use energyFilter_class,   only : energyFilter
  use testFilter_class,     only : testFilter


  implicit none
  private

  public :: new_tallyFilter

  ! *** ADD NAME OF A NEW TALLY FILTER HERE ***!
  ! List that contains all accaptable types of tallyFilters
  ! It is printed if type was unrecognised
  ! NOTE:
  ! For now  it is necessary to adjust trailing blanks so all enteries have the same length
  character(nameLen),dimension(*),parameter :: AVALIBLE_tallyFilters = [ 'energyFilter']

contains

  !!
  !! Allocate new allocatable tallyFilter to a specific type
  !! If new is allocated it deallocates it
  !!
  subroutine new_tallyFilter(new,dict)
    class(tallyFilter),allocatable, intent(inout) :: new
    class(dictionary), intent(in)                 :: dict
    character(nameLen)            :: type
    character(100),parameter      :: Here = 'new_tallyFilter (tallyFilterFactory_func.f90)'

    ! Deallocate new if allocated
    if(allocated(new)) deallocate(new)

    ! Obtain string that specifies type to be built
    call dict % get(type,'type')

    ! Allocate approperiate subclass of tallyFilter
    ! *** ADD CASE STATEMENT FOR A NEW TALLY MAP BELOW ***!
    select case(type)
      case('energyFilter')
        allocate(energyFilter :: new)
        call new % init(dict)

      case('testFilter')
        allocate(testFilter :: new)
        call new % init(dict)

     !*** NEW TALLY MAP TEMPLATE ***!
     !case('<newtallyFilterName>')
     !  allocate(<newtallyFilterName> :: new)
     !  call new % init(dict)
     !
      case default
        print *, AVALIBLE_tallyFilters
        call fatalError(Here, 'Unrecognised type of tallyFilter: ' // trim(type))

    end select

  end subroutine new_tallyFilter

end module tallyFilterFactory_func
