module collisionProcessorFactory_func

  use numPrecision
  use genericProcedures, only : fatalError
  use dictionary_class,  only : dictionary

  ! Abstract interface
  use collisionProcessor_inter, only : collisionProcessor

  ! Implementation
  use neutronCEstd_class, only : neutronCEstd
  use neutronCEimp_class, only : neutronCEimp
  use neutronMGstd_class, only : neutronMGstd

  implicit none
  private

  public :: new_collisionProcessor

  ! *** ADD NAME OF A NEW COLLISION PROCESSOR HERE ***!
  ! List that contains all accaptable types of collisionProcessors
  ! It is printed if type was unrecognised
  ! NOTE:
  ! For now  it is necessary to adjust trailing blanks so all enteries have the same length
  character(nameLen),dimension(*),parameter :: AVALIBLE_collisionProcessors = [ 'neutronCEstd',&
                                                                                'neutronCEimp',&
                                                                                'neutronMGstd']

contains

  !!
  !! Allocate new allocatable collisionProcessor to a specific type
  !! If new is allocated it deallocates it
  !!
  subroutine new_collisionProcessor(new,dict)
    class(collisionProcessor),allocatable, intent(inout) :: new
    class(dictionary), intent(in)                        :: dict
    character(nameLen)                                   :: type
    character(100),parameter      :: Here = 'new_collisionProcessor (collisionProcessorFactory_func.f90)'

    ! Deallocate new if allocated
    if(allocated(new)) deallocate(new)

    ! Obtain string that specifies type to be built
    call dict % get(type,'type')

    ! Allocate approperiate subclass of collisionProcessor
    ! *** ADD CASE STATEMENT FOR A NEW COLLISION PROCESSOR BELOW ***!
    select case(type)
      case('neutronCEstd')
        allocate(neutronCEstd :: new)
        call new % init(dict)

      case('neutronCEimp')
        allocate(neutronCEimp :: new)
        call new % init(dict)

      case('neutronMGstd')
        allocate(neutronMGstd :: new)
        call new % init(dict)

     !*** NEW COLLISION PROCESSOR TEMPLATE ***!
     !case('<newcollisionProcessorName>')
     !  allocate(<newcollisionProcessorName> :: new)
     !  call new % init(dict)
     !
      case default
        print *, AVALIBLE_collisionProcessors
        call fatalError(Here, 'Unrecognised type of collisionProcessor: ' // trim(type))

    end select

  end subroutine new_collisionProcessor

end module collisionProcessorFactory_func
