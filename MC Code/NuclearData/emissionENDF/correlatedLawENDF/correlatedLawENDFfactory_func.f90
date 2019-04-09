module correlatedLawENDFfactory_func

  use numPrecision
  use endfConstants
  use genericProcedures,       only : fatalError, numToChar
  use aceCard_class,           only : aceCard

  ! Correlated Laws
  use correlatedLawENDF_inter, only : correlatedLawENDF
  use kalbach87_class,         only : kalbach87
  use endfLaw61_class,         only : endfLaw61
  use nBodyPhaseSpace_class,   only : nBodyPhaseSpace

  implicit none
  private

  public :: new_correlatedLawENDF

contains

  !!
  !! Allocates a new  correlatedLawENDF from aceCard and MT number
  !! aceCard can be in any poistion. Its position changes at output.
  !!
  subroutine new_correlatedLawENDF(new, ACE, MT)
    class(correlatedLawENDF),allocatable, intent(inout) :: new
    type(aceCard),intent(inout)                         :: ACE
    integer(shortInt), intent(in)                       :: MT
    integer(shortInt)                                   :: LOCB, LNW, LAW, loc
    real(defReal)                                       :: Q, A
    character(100),parameter :: Here='new_correlatedLawENDF (correlatedLawENDFfactory_func.f90)'

    ! Deallocate new if allocated
    if(allocated(new)) deallocate(new)

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
      call fatalError(Here,'Multiple correlated laws for nuclide:' // trim(ACE % ZAID) //' MT:' &
                            // numToChar(MT))

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

      case(endfEnergyLaw61)
        allocate(new, source = endfLaw61(ACE))

      case(nBodyPhaseSpaceDistribution)
        !! Get Q & A value
        Q = ACE % QforMT(MT)
        A = ACE % AW
        allocate(new, source = nBodyPhaseSpace(ACE, Q, A))

      case default
        print *, 'Energy Law Type :', LAW
        call fatalError(Here,'Correlated Law Type is not recognised or yet supported ')

    end select

  end subroutine new_correlatedLawENDF

end module correlatedLawENDFfactory_func
