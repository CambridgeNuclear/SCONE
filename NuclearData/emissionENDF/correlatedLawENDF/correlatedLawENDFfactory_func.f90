module correlatedLawENDFfactory_func

  use numPrecision
  use endfConstants
  use genericProcedures,       only : fatalError, numToChar
  use aceCard_class,           only : aceCard

  ! Correlated Laws
  use correlatedLawENDF_inter,      only : correlatedLawENDF
  use kalbach87_class,              only : kalbach87
  use endfLaw61_class,              only : endfLaw61
  use nBodyPhaseSpace_class,        only : nBodyPhaseSpace
  use multipleCorrelatedLaws_class, only : multipleCorrelatedLaws

  implicit none
  private

  public :: new_correlatedLawENDF

contains

  !!
  !! Allocate a new correlated angle-energy law from ACE card
  !!
  !! Args:
  !!   ACE [inout] -> ACE card with data. Can be set to any place.
  !!   MT [in]     -> MT number of the requested reaction
  !!
  !! Errors:
  !!   FatalError if MT reaction is not correlated (LOCB /= -1)
  !!
  subroutine new_correlatedLawENDF(new, ACE, MT)
    class(correlatedLawENDF),allocatable, intent(inout) :: new
    type(aceCard),intent(inout)                         :: ACE
    integer(shortInt), intent(in)                       :: MT
    integer(shortInt)                                   :: LOCB, LNW, LAW, loc, root
    integer(shortInt)                                   :: N, i, NR, NEne
    class(multipleCorrelatedLaws),allocatable           :: multiLaw
    character(100),parameter :: Here = 'new_correlatedLawENDF (correlatedLawENDFfactory_func.f90)'

    ! Deallocate new if allocated
    if(allocated(new)) deallocate(new)

    ! Verify that the energy law is indeed coreelated
    LOCB = ACE % LOCBforMT(MT)

    if(LOCB /= LOCB_CORRELATED) then
      call fatalError(Here,'Reaction under MTdoes not have correlated mu-energy distribution')
    end if

    ! Set aceCard read head to beginning of energy data for MT
    call ACE % setToEnergyMT(MT)
    root = ACE % getRootAddress('energyLawsMT')

    ! Read location of next energy law. If LNW == 0 only one law is given
    LNW = ACE % readInt()

    if (LNW == 0) then
      LAW = ACE % readInt()
      loc = ACE % readInt()
      call buildENDFLaw(new, LAW, MT, root, loc, ACE)

    else

      ! Calculate number of correlated laws
      N = 1
      do while (LNW /= 0)
        N = N + 1
        call ACE % setRelativeTo(root, LNW)
        LNW = ACE % readINT()
        if (N > 100) call fatalError(Here, 'Infinate loop in LAW count. Terminating.')
      end do

      ! Allocate space
      allocate(multiLaw)
      call multiLaw % init(N)

      ! Reset to the begining
      call ACE % setToEnergyMT(MT)
      LNW = ACE % readInt()

      ! Read all laws
      do i = 1, N
        LAW = ACE % readInt()
        loc = ACE % readInt()

        ! Read number of interpolation regions
        NR = ACE % readInt()
        associate(bounds => ACE % readIntArray(NR), interENDF => ACE % readIntArray(NR))
          ! Read number of energy regions
          NEne = ACE % readInt()
          associate( eGrid => ACE % readRealArray(NEne), pdf => ACE % readRealArray(NEne))
            ! Build energy law
            call buildENDFLaw(new, LAW, MT, root, loc, ACE)

            if(NR == 0) then
              call multiLaw % addLaw(new, eGrid, pdf)
            else
              call multiLaw % addLaw(new, eGrid, pdf, bounds, interENDF)
            end if
          end associate
        end associate

        ! Move to the next low and read new LNW
        ! Protect against last iteration
        if (LNW /= 0) then
          call ACE % setRelativeTo(root, LNW)
          LNW = ACE % readInt()
        end if
      end do

      ! Verify that all laws were read
      if(LNW /= 0) call fatalError(Here,'LNW is not 0 after reading all correlated laws. It is ' // &
                                         numToChar(LNW) //' Somthing failed')
      ! Move finished multiple laws to new
      call move_alloc(multiLaw, new)

    end if

  end subroutine new_correlatedLawENDF

  !!
  !! Helper procedure to alocate a correlated energy law
  !!
  !! See MCNP 4 Manual Appendix F Table F-14 for extra info
  !!
  !! Args:
  !!   lawENDF [inout] -> polymorphic correlated law to be allocated
  !!   LAW [in]        -> Integer ID of the correlated law
  !!   MT [in]         -> MT number of the reaction
  !!   root [in]       -> Root location on XSS in ACE for this correlated law
  !!   offset [in]     -> Location relative to root on XSS in ACE for the law
  !!   ACE [inout]     -> ACE Card with data
  !!
  !! Errors:
  !!   Will crash if root and offset point to incorrect location
  !!
  subroutine buildENDFLaw(lawENDF, LAW, MT, root, offset, ACE)
    class(correlatedLawENDF), allocatable, intent(inout) :: lawENDF
    integer(shortInt), intent(in)                        :: LAW
    integer(shortInt), intent(in)                        :: MT
    integer(shortInt), intent(in)                        :: root
    integer(shortInt), intent(in)                        :: offset
    type(aceCard), intent(inout)                         :: ACE
    real(defReal)                                        :: Q, A
    character(100),parameter :: Here = 'buildENDFLaw (correlatedLawENDFfactory_func.f90)'

    ! Deallocate lawENDF if allocated
    if(allocated(lawENDF)) deallocate(lawENDF)

    ! Set ACE card to correct location and build the data
    call ACE % setRelativeTo(root, offset)

    select case(LAW)
      case(kalbach87Formalism)
        allocate(lawENDF, source = kalbach87(ACE))

      case(endfEnergyLaw61)
        allocate(lawENDF, source = endfLaw61(ACE))

      case(nBodyPhaseSpaceDistribution)
        !! Get Q & A value
        Q = ACE % QforMT(MT)
        A = ACE % AW
        allocate(lawENDF, source = nBodyPhaseSpace(ACE, Q, A))

      case default
        call fatalError(Here,'Correlated Law Type is not recognised or yet &
                              &supported: '//numToChar(LAW))

    end select

  end subroutine buildENDFLaw

end module correlatedLawENDFfactory_func
