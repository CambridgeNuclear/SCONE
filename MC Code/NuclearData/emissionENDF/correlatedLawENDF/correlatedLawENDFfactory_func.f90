module correlatedLawENDFfactory_func

  use numPrecision
  use endfConstants
  use genericProcedures,       only : fatalError
  use aceCard_class,           only : aceCard

  ! Correlated Laws
  use correlatedLawENDF_inter, only : correlatedLawENDF
  use kalbach87_class,         only : kalbach87

  implicit none
  private

  public :: new_correlatedLawENDF

contains

  !!
  !! Returns an allocatable correlatedLawENDF from aceCard and MT number
  !! aceCard can be in any poistion. Its position changes at output.
  !!
  function new_correlatedLawENDF(ACE,MT) result(new)
    type(aceCard),intent(inout)          :: ACE
    integer(shortInt), intent(in)        :: MT
    class(correlatedLawENDF),allocatable :: new
    integer(shortInt)                    :: LOCB,LNW,LAW,loc
    character(100),parameter :: Here='new_correlatedLawENDF (correlatedLawENDFfactory_func.f90)'

    ! Verify that the energy law is indeed coreelated
    LOCB = ACE % LOCBforMT(MT)

    if(LOCB /= LOCB_CORRELATED) then
      call fatalError(Here,'Reaction under MTdoes not have correlated mu-energy distribution')
    end if

    ! Set aceCard read head to beginning of energy data for MT
    call ACE % setToEnergyMT(MT)

    ! Read location of next energy law. If LNW == 0 only one law is given
    LNW = ACE % readInt()

    ! Give error if multiple laws are present
    if (LNW /= 0) then
      call fatalError(Here,'Multiple energy laws for a single MT are not yet supported')

    end if

    ! Read energy Law type and location
    LAW = ACE % readInt()
    loc = ACE % readInt()

    ! Set ACE to location of the LAW
    call ACE % setToEnergyLaw(loc)

    ! Allocate new object
    select case(LAW)
      case(kalbach87Formalism)
        allocate(new, source = kalbach87(ACE))

      case default
        print *, 'Energy Law Type :', LAW
        call fatalError(Here,'Correlated Law Type is not recognised or yet supported ')

    end select

  end function new_correlatedLawENDF

  !!
  !! Returns a pointer to allocated correlatedLawENDF from aceCard and MT number
  !!
  function new_correlatedLawENDF_ptr(ACE,MT) result(new)
    type(aceCard), intent(inout)     :: ACE
    integer(shortInt), intent(in)    :: MT
    class(correlatedLawENDF),pointer :: new

    ! Allocate pointer and copy data from local allocatable
    allocate(new, source = new_correlatedLawENDF(ACE,MT))

  end function new_correlatedLawENDF_ptr

end module correlatedLawENDFfactory_func
