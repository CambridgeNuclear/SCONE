module asciiOutputFactory_func

  use numPrecision
  use genericProcedures, only : fatalError

  ! Interface
  use asciiOutput_inter, only : asciiOutput

  ! Implementations
  use asciiMATLAB_class,  only : asciiMATLAB
  use dummyPrinter_class, only : dummyPrinter

  implicit none
  private

  ! *** ADD NAME OF A NEW ASCII OUTPUT HERE ***!
  ! List that contains all accaptable types of ascii Output printers
  ! It is printed if type was unrecognised
  ! NOTE:
  ! For now  it is necessary to adjust trailing blanks so all enteries have the same length
  character(nameLen),dimension(*),parameter :: AVALIBLE_asciiOutputs = [ 'asciiMATLAB ',&
                                                                         'dummyPrinter']

  public :: new_asciiOutput

contains

  !!
  !! Returns an instance of asciiOutput from type character
  !!
  function new_asciiOutput(type) result(new)
    character(nameLen), intent(in) :: type
    class(asciiOutput),allocatable :: new
    character(100),parameter :: Here = 'new_asciiOutput (asciiOutputFactory_func.f90)'

    ! Allocate approperiate subclass of asciiOutput
    ! *** ADD CASE STATEMENT FOR A NEW ASCII OUTPUT BELOW ***!
    select case(type)
      case('asciiMATLAB')
        allocate(new, source = asciiMATLAB() )

      case('dummyPrinter')
        allocate(new, source = dummyPrinter() )

      case default
        print *, AVALIBLE_asciiOutputs
        call fatalError(Here, 'Unrecognised type of asciiOutput: ' // trim(type))

    end select

  end function new_asciiOutput

end module asciiOutputFactory_func
