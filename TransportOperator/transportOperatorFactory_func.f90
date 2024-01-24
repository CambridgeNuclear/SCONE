!!
!! Module to build transport operators
!!
module transportOperatorFactory_func

  use numPrecision
  use genericProcedures, only : fatalError
  use dictionary_class,  only : dictionary

  ! Transport Operators
  use transportOperator_inter,          only : transportOperator
  use transportOperatorST_class,        only : transportOperatorST
  use transportOperatorDT_class,        only : transportOperatorDT
  use transportOperatorHT_class,        only : transportOperatorHT
  !use transportOperatorDynamicDT_class, only : transportOperatorDynamicDT

  implicit none
  private

  ! List that contains all accaptable types of transport operators
  ! It is printed if type was unrecognised
  ! NOTE:
  ! For now  it is necessary to adjust trailing blanks so all enteries have the same length
  character(nameLen),dimension(*),parameter :: AVALIBLE_transportOps = [ 'transportOperatorST', &
                                                                         'transportOperatorDT', &
                                                                         'transportOperatorHT']

  public :: new_transportOperator

contains

  !!
  !! Allocate new allocatable transportOperator to a specific type
  !! If new is allocated it deallocates it
  !!
  subroutine new_transportOperator(new, dict)
    class(transportOperator),allocatable, intent(inout):: new
    class(dictionary), intent(in)                      :: dict
    character(nameLen)                                 :: type
    character(100),parameter :: Here = 'new_transportOperator (transportOperatorFactory_func.f90)'

    if(allocated(new)) deallocate(new)

    ! Obtain string that specifies type to be built
    call dict % get(type,'type')

    ! Allocate approperiate subclass of transportOperator
    select case(type)
      case('transportOperatorST')
        allocate( transportOperatorST :: new)

      case('transportOperatorDT')
        allocate( transportOperatorDT :: new)

      case('transportOperatorHT')
        allocate( transportOperatorHT :: new)

      case default
        print *, AVALIBLE_transportOps
        call fatalError(Here, 'Unrecognised type of transportOperator: ' // trim(type))

    end select

    ! Initialise new transport operator
    call new % init(dict)

  end subroutine new_transportOperator


end module transportOperatorFactory_func
