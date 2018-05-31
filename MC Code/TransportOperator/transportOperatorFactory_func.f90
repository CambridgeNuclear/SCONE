!!
!! Module to build transport operators
!!
module transportOperatorFactory_func

  use numPrecision
  use genericProcedures, only : fatalError
  use dictionary_class,  only : dictionary

  ! Nuclear Data an geometry interface
  use nuclearData_inter, only : nuclearData
  use geometry_class,    only : geometry

  ! Transport Operators
  use transportOperator_inter,   only : transportOperator
  use transportOperatorST_class, only : transportOperatorST
  use transportOperatorDT_class, only : transportOperatorDT

  implicit none
  private

  ! *** ADD NAME OF A NEW TRANSPORT OPERATOR HERE ***!
  ! List that contains all accaptable types of transport operators
  ! It is printed if type was unrecognised
  ! NOTE:
  ! For now  it is necessary to adjust trailing blanks so all enteries have the same length
  character(nameLen),dimension(*),parameter :: AVALIBLE_transportOps = [ 'transportOperatorST', &
                                                                         'transportOperatorDT']

  public :: new_transportOperator_ptr

contains

  function new_transportOperator_ptr(nucData,geom,dict) result(new)
    class(nuclearData),pointer,intent(in) :: nucData
    class(geometry),pointer,intent(in)    :: geom
    class(dictionary), intent(in)     :: dict
    class(transportOperator),pointer  :: new
    character(nameLen)                :: type
    character(100),parameter :: Here = 'new_transportOperator_ptr (transportOperatorFactory_func.f90)'

    ! Obtain string that specifies type to be built
    call dict % get(type,'type')

    ! Allocate approperiate subclass of transportOperator
    ! *** ADD CASE STATEMENT FOR A NEW TRANSPORT OPERATOR BELOW ***!
    select case(type)
      case('transportOperatorST')
        ! Allocate and initialise
        allocate( transportOperatorST :: new)
        call new % init(nucData, geom, dict)

      case('transportOperatorDT')
        ! Allocate and initialise
        allocate( transportOperatorDT :: new)
        call new % init(nucData, geom, dict)

      case default
        print *, AVALIBLE_transportOps
        call fatalError(Here, 'Unrecognised type of transportOperator: ' // trim(type))

    end select


  end function new_transportOperator_ptr


end module transportOperatorFactory_func
