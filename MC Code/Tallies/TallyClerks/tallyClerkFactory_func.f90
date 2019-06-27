module tallyClerkFactory_func

  use numPrecision
  use genericProcedures, only : fatalError
  use dictionary_class,  only : dictionary

  ! tallyClerk interface
  use tallyClerk_inter,    only : tallyClerk

  ! tallyClerk implementations
  use keffAnalogClerk_class,     only : keffAnalogClerk
  use keffImplicitClerk_class,   only : keffImplicitClerk
  use collisionClerk_class,      only : collisionClerk
  use simpleFMClerk_class,       only : simpleFMClerk
  use dancoffBellClerk_class,    only : dancoffBellClerk
  use shannonEntropyClerk_class, only : shannonEntropyClerk
  use centreOfMassClerk_class,   only : centreOfMassClerk

  implicit none
  private

  public :: new_tallyClerk

  ! *** ADD NAME OF A NEW TALLY FILTER HERE ***!
  ! List that contains all accaptable types of tallyClerks
  ! It is printed if type was unrecognised
  ! NOTE:
  ! For now  it is necessary to adjust trailing blanks so all enteries have the same length
  character(nameLen),dimension(*),parameter :: AVALIBLE_tallyClerks = [ 'keffAnalogClerk    ',&
                                                                        'keffImplicitClerk  ',&
                                                                        'collisionClerk     ',&
                                                                        'simpleFMClerk      ',&
                                                                        'dancoffBellClerk   ',&
                                                                        'shannonEntropyClerk',&
                                                                        'centreOfMassClerk  ']

contains

  !!
  !! Allocate new allocatable tallyClerk to a specific type
  !! If new is allocated it deallocates it
  !!
  subroutine new_tallyClerk(new, dict, name)
    class(tallyClerk),allocatable, intent(inout) :: new
    class(dictionary), intent(in)                :: dict
    character(nameLen),intent(in)                :: name
    character(nameLen)            :: type
    character(100),parameter      :: Here = 'new_tallyClerk (tallyClerkFactory_func.f90)'

    ! Deallocate new if allocated
    if(allocated(new)) deallocate(new)

    ! Obtain string that specifies type to be built
    call dict % get(type,'type')

    ! Allocate approperiate subclass of tallyClerk
    ! *** ADD CASE STATEMENT FOR A NEW TALLY MAP BELOW ***!
    select case(type)
     case('keffAnalogClerk')
       allocate(keffAnalogClerk :: new)
       call new % init(dict, name)

     case('keffImplicitClerk')
       allocate(keffImplicitClerk :: new)
       call new % init(dict, name)

     case('collisionClerk')
       allocate(collisionClerk :: new)
       call new % init(dict, name)

     case('simpleFMClerk')
       allocate(simpleFMClerk :: new)
       call new % init(dict, name)

     case('dancoffBellClerk')
       allocate(dancoffBellClerk :: new)
       call new % init(dict, name)

     case('shannonEntropyClerk')
       allocate(shannonEntropyClerk :: new)
       call new % init(dict, name)

     case('centreOfMassClerk')
       allocate(centreOfMassClerk :: new)
       call new % init(dict, name)

     !*** NEW TALLY MAP TEMPLATE ***!
     !case('<newtallyClerkName>')
     !  allocate(<newtallyClerkName> :: new)
     !  call new % init(dict, name)
     !
      case default
        print *, AVALIBLE_tallyClerks
        call fatalError(Here, 'Unrecognised type of tallyClerk: ' // trim(type))

    end select

  end subroutine new_tallyClerk

end module tallyClerkFactory_func
