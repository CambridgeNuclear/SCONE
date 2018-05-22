module tallyClerkFactory_func

  use numPrecision
  use dictionary_class,        only : dictionary
  use tallyClerk_inter,        only : tallyClerk
  use keffClerk_inter,         only : keffClerk
  use keffActiveClerk_class,   only : keffActiveClerk
  use keffInactiveClerk_class, only : keffInactiveClerk

  implicit none
  private

  interface new_tallyClerk
    module procedure newTallyCLerc_alloc
  end interface

  ! *** ADD NAME OF A NEW TALLY CLERK HERE ***!
  ! List that contains all accaptable types of tally Clerks
  ! It is printed if type was unrecognised
  character(nameLen),dimension(:),parameter :: avalible_tallyClerks = [ 'keffActiveClerk',&
                                                                        'keffInactiveClerk' ]

contains

  subroutine new_tallyClerk_alloc(tallyClerk,dict)
    class(tallyClerk),allocatable,intent(inout) :: self
    class(dictionary), intent(in)               :: dict
    character(nameLen)                          :: type

    ! Obtain string that specifies type to be built
    call dict % get(type,'type')

    select case(type)
      case('keffActiveClerk')

      case('keffInactiveClerk')

      case default

    end select


  end subroutine newTallyClerk_alloc



end module tallyClerkFactory_func
