!!
!! Module to build pieceConstantField
!!
module pieceConstantFieldFactory_func

  use numPrecision
  use errors_mod,          only : fatalError
  use dictionary_class,    only : dictionary

  use pieceConstantField_inter, only : pieceConstantField
  use cartesianField_class,     only : cartesianField

  implicit none
  private

  ! List that contains all acceptable types of pieceConstantFields
  ! It is printed if type was unrecognised
  ! NOTE:
  ! For now  it is necessary to adjust trailing blanks so all enteries have the same length
  character(nameLen),dimension(*),parameter :: AVAILABLE_fields = [ 'cartesianField']

  public :: new_pieceConstantField

contains

  !!
  !! Allocate new allocatable pieceConstantField to a specific type
  !! If new is allocated it deallocates it
  !!
  subroutine new_pieceConstantField(new, dict)
    class(pieceConstantField),allocatable, intent(inout):: new
    class(dictionary), intent(in)                       :: dict
    character(nameLen)                                  :: type
    character(100),parameter :: Here = 'new_pieceConstantField (pieceConstantFieldFactory_func.f90)'

    if(allocated(new)) deallocate(new)

    ! Obtain string that specifies type to be built
    call dict % get(type,'type')

    ! Allocate approperiate subclass of pieceConstantField
    select case(type)
      case('cartesianField')
        allocate( cartesianField :: new)

      case default
        print *, AVAILABLE_fields
        call fatalError(Here, 'Unrecognised type pieceConstantField: ' // trim(type))

    end select

    ! Initialise new pieceConstantField
    call new % init(dict)

  end subroutine new_pieceConstantField


end module pieceConstantFieldFactory_func
