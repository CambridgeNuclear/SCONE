module tallyClerkFactory_func

  use numPrecision
  use genericProcedures, only : fatalError
  use dictionary_class,  only : dictionary

  ! tallyClerk interface
  use tallyClerk_inter,    only : tallyClerk

  ! tallyClerk implementations
  use keffAnalogClerk_class,           only : keffAnalogClerk
  use keffImplicitClerk_class,         only : keffImplicitClerk
  use collisionClerk_class,            only : collisionClerk
  use collisionProbabilityClerk_class, only : collisionProbabilityClerk
  use trackClerk_class,                only : trackClerk
  use simpleFMClerk_class,             only : simpleFMClerk
  use dancoffBellClerk_class,          only : dancoffBellClerk
  use shannonEntropyClerk_class,       only : shannonEntropyClerk
  use centreOfMassClerk_class,         only : centreOfMassClerk
  use energyWeightClerk_class,         only : energyWeightClerk
  use mgXsClerk_class,                 only : mgXsClerk

  implicit none
  private

  public :: new_tallyClerk

  ! List that contains all accaptable types of tallyClerks
  ! It is printed if type was unrecognised
  ! NOTE:
  ! For now  it is necessary to adjust trailing blanks so all enteries have the same length
  character(nameLen),dimension(*),parameter :: AVALIBLE_tallyClerks = [ 'keffAnalogClerk          ',&
                                                                        'keffImplicitClerk        ',&
                                                                        'collisionClerk           ',&
                                                                        'collisionProbabilityClerk',&
                                                                        'trackClerk               ',&
                                                                        'simpleFMClerk            ',&
                                                                        'shannonEntropyClerk      ',&
                                                                        'centreOfMassClerk        ',&
                                                                        'dancoffBellClerk         ',&
                                                                        'energyWeightClerk        ']
                                                                        'mgXsClerk                ']

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
    select case(type)
     case('keffAnalogClerk')
       allocate(keffAnalogClerk :: new)

     case('keffImplicitClerk')
       allocate(keffImplicitClerk :: new)

     case('collisionClerk')
       allocate(collisionClerk :: new)

     case('collisionProbabilityClerk')
       allocate(collisionProbabilityClerk :: new)

     case('trackClerk')
       allocate(trackClerk :: new)

     case('simpleFMClerk')
       allocate(simpleFMClerk :: new)

     case('dancoffBellClerk')
       allocate(dancoffBellClerk :: new)

     case('shannonEntropyClerk')
       allocate(shannonEntropyClerk :: new)

     case('centreOfMassClerk')
       allocate(centreOfMassClerk :: new)

     case('energyWeightClerk')
       allocate(energyWeightClerk :: new)

     case('mgXsClerk')
       allocate(mgXsClerk :: new)

     case default
        print *, AVALIBLE_tallyClerks
        call fatalError(Here, 'Unrecognised type of tallyClerk: ' // trim(type))

    end select

    ! Initialise new clerk
    call new % init(dict, name)

  end subroutine new_tallyClerk

end module tallyClerkFactory_func
