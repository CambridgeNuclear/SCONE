module tallyClerkFactory_func

  use numPrecision
  use genericProcedures,       only : fatalError
  use dictionary_class,        only : dictionary

  ! Tally Clerks
  use tallyClerk_inter,        only : tallyClerk
  use keffActiveClerk_class,   only : keffActiveClerk
  use keffInactiveClerk_class, only : keffInactiveClerk

  implicit none
  private

  public :: new_tallyClerk
  public :: new_tallyClerk_ptr

  ! *** ADD NAME OF A NEW TALLY CLERK HERE ***!
  ! List that contains all accaptable types of tally Clerks
  ! It is printed if type was unrecognised
  ! NOTE:
  ! For now  it is necessary to adjust trailing blanks so all enteries have the same length
  character(nameLen),dimension(*),parameter :: AVALIBLE_tallyClerks = [ 'keffActiveClerk  ',&
                                                                        'keffInactiveClerk' ]


contains

  ! *** ADD ENTERY TO THIS PROCEDURE FOR A NEW TALLY CLERK ***
  !!
  !! Allocates an allocatable tallyClerk based on dictionary
  !!
  function new_tallyClerk(dict,name) result(new)
    class(dictionary), intent(in)               :: dict
    character(nameLen), intent(in)              :: name
    class(tallyClerk),allocatable               :: new
    character(nameLen)                          :: type
    character(100),parameter :: Here = 'new_tallyClerk (tallyClerkFactory_func.f90)'

    ! Obtain string that specifies type to be built
    call dict % get(type,'type')

    ! Allocate approperiate subclass of tallyClerk
    ! *** ADD CASE STATEMENT FOR A NEW TALLY CLERK BELOW ***!
    select case(type)
      case('keffActiveClerk')
        allocate(new, source = keffActiveClerk(dict,name) )

      case('keffInactiveClerk')
        allocate(new, source = keffInactiveClerk(dict,name) )

     !*** NEW TALLY CLERK TEMPLATE ***!
     !case('<newTallyClerkName>')
     !  allocate(new, source = <newTallyClerkName>(dict) )
     !
      case default
        print *, AVALIBLE_tallyClerks
        call fatalError(Here, 'Unrecognised type of tallyClerk: ' // trim(type))

    end select

  end function new_TallyClerk

  !!
  !! Allocates a null pointer to tallyClerk to an instance based on dictionary
  !! If pointer is alrady allocated it returns an error
  !!
  function new_tallyClerk_ptr(dict,name) result(new)
    class(dictionary), intent(in)            :: dict
    character(nameLen), intent(in)           :: name
    class(tallyClerk),pointer                :: new
    character(100),parameter :: Here = 'new_tallyClerk_ptr (tallyClerkFactory_func.f90)'

    ! Allocate pointer and copy data from local allocatable
    allocate( new, source = new_TallyClerk(dict, name) )

  end function new_tallyClerk_ptr

end module tallyClerkFactory_func
