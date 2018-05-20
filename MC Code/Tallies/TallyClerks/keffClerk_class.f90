module keffClerk_class

  use numPrecision
  use tallyCodes
  use genericProcedures,          only : fatalError
  use particle_class,             only : particle, phaseCoord
  use particleDungeon_class,      only : particleDungeon
  use tallyClerk_inter,           only : tallyClerk
  use tallyEstimator_class,       only : tallyScore, tallyCounter
  use transportNuclearData_inter, only : transportNuclearData

  use xsMacroSet_class,           only : xsMacroSet_ptr

  implicit none
  private

  type, public,extends(tallyClerk) :: keffClerk
    private

    integer(shortInt)                    :: cycleCount = 0

    type(tallyScore)                     :: impProd     ! Implicit neutron production
    type(tallyScore)                     :: impAbs      ! Implicit neutron absorbtion
    type(tallyScore)                     :: anaLeak     ! Analog neutron leakage

    type(tallyCounter)                   :: k_analog
    type(tallyCounter)                   :: k_imp

    real(defReal)                        :: startWgt

  contains
    ! Deferred Interface Procedures
    procedure :: validReports
    procedure :: display

    ! Overwrite report procedures
    procedure :: reportInColl
    procedure :: reportHist
    procedure :: reportCycleStart
    procedure :: reportCycleEnd
  end type keffClerk

contains

  !!
  !! Return codes for reports this clerk accepts
  !!
  function validReports(self) result(validCodes)
    class(keffClerk), intent(in)               :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    validCodes = [inColl_CODE, cycleStart_CODE, cycleEnd_CODE, hist_CODE]

  end function validReports

  !!
  !! Display progress current estimate of k-eff with STD
  !!
  subroutine display(self)
    class(keffClerk), intent(in) :: self
    real(defReal)                :: k_imp, k_analog, STD_imp, STD_analog

    ! Obtain current estimates of k analog and implicit
    call self % k_imp % getEstimate(k_imp, STD_imp, self % cycleCount)
    call self % k_analog % getEstimate(k_analog, STD_analog, self % cycleCount)

    ! Print estimates to a console
    print '(A,F8.5,A,F8.5)', 'k-eff (implicit): ', k_imp, ' +/- ', STD_imp
    print '(A,F8.5,A,F8.5)', 'k-eff (analog): ',  k_analog, ' +/- ', STD_analog

  end subroutine display

  !!
  !! Process collision report
  !!
  subroutine reportInColl(self,p)
    class(keffClerk), intent(inout)       :: self
    class(particle), intent(in)           :: p
    type(xsMacroSet_ptr)                  :: XSs
    real(defReal)                         :: totalXS, nuFissXS, absXS, flux
    real(defReal)                         :: s1, s2
    character(100), parameter  :: Here = 'reportInColl (keffClerk_class.f90)'

    ! Obtain XSs
    ! Check if it dynamic type is supported
    ! If it is obtain macroscopic XSs
    ! It it isn't throw error
    associate (xsData => p % xsData)
      select type(xsData)
        class is (transportNuclearData)
          call xsData % getMatMacroXS(XSs, p, p % matIdx)

        class default
          call fatalError(Here,'Dynamic type of XS data attached to particle is not transportNuclearData')

      end select
    end associate

    totalXS  = XSs % totalXS()
    nuFissXS = XSs % nuFissionXS()
    absXS    = XSs % captureXS() + XSs % fissionXS()

    ! Calculate flux and scores
    flux = p % w / totalXS

    s1 = nuFissXS * flux
    s2 = absXS * flux

    ! Add scores to counters
    call self % impProd % add(s1)
    call self % impAbs  % add(s2)

  end subroutine reportInColl

  !!
  !! Process history report
  !!
  subroutine reportHist(self,pre,post,fate)
    class(keffClerk), intent(inout)      :: self
    class(phaseCoord), intent(in)        :: pre
    class(particle), intent(in)          :: post
    integer(shortInt),intent(in)         :: fate
    real(defReal)   :: histWgt


    if( fate == leak_FATE) then
      ! Obtain and score history weight
      histWgt = pre % wgt

      ! Score analog leakage
      call self % anaLeak % add(histWgt)

    end if

  end subroutine reportHist

  !!
  !! Process beginning of a cycle
  !!
  subroutine reportCycleStart(self,start)
    class(keffClerk), intent(inout)      :: self
    class(particleDungeon), intent(in)   :: start

    self % startWgt  = start % popWeight()

  end subroutine reportCycleStart

  !!
  !! Process end of the cycle
  !!
  subroutine reportCycleEnd(self,end)
    class(keffClerk), intent(inout)      :: self
    class(particleDungeon), intent(in)   :: end
    real(defReal)                        :: endWgt, k_est
    real(defReal)      :: nuFiss, absorb, leakage, collCount, histCount, k_cycle

    ! Obtain end of cycle weight and k value used to change fission site generation rate
    endWgt = end % popWeight()
    k_cycle = end % k_eff

    ! Calculate and score analog estimate of k-eff
    k_est =  endWgt / self % startWgt * k_cycle
    call self % k_analog % addEstimate(k_est)

    ! Calculate and score implicit estimate of k_eff
    nuFiss  = self % impProd % get()
    absorb  = self % impAbs % get()
    leakage = self % anaLeak % get()

    k_est = nuFiss / (absorb + leakage  )

    call self % k_imp % addEstimate(k_est)

    ! Reset score counters
    call self % impProd % reset()
    call self % impAbs % reset()
    call self % anaLeak % reset()

    ! Increas counter of cycles
    self % cycleCount = self % cycleCount + 1

  end subroutine reportCycleEnd

end module keffClerk_class
