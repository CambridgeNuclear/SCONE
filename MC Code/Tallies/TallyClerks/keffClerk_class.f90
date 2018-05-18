module keffClerk_class

  use numPrecision
  use genericProcedures,          only : fatalError
  use particle_class,             only : particle, phaseCoord
  use particleDungeon_class,      only : particleDungeon
  use tallyClerk_inter,           only : tallyClerk,  inColl_CODE, outColl_CODE, path_CODE, &
                                         trans_CODE, hist_CODE, cycleStart_CODE, cycleEnd_CODE
  use tallyCounter_class,         only : tallyCounter
  use transportNuclearData_inter, only : transportNuclearData

  use xsMacroSet_class,           only : xsMacroSet_ptr

  implicit none
  private

  type, public,extends(tallyClerk) :: keffClerk
    private

    integer(shortInt)                    :: cycleCount = 0       !
    type(tallyCounter)                   :: collCount   ! Total collision weight
    type(tallyCounter)                   :: histCount        ! Total histories weight


    type(tallyCounter)                   :: impProd     ! Implicit neutron production
    type(tallyCounter)                   :: impAbs      ! Implicit neutron absorbtion
    type(tallyCounter)                   :: anaLeak     ! Analog neutron leakage

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
    ! Implement display

    call self % k_imp % getScore(k_imp, STD_imp, self % cycleCount)
    call self % k_analog % getScore(k_analog, STD_analog, self % cycleCount)

    print *, 'k-eff (implicit): ', k_imp, ' +/- ', STD_imp
    print *, 'k-eff (analog): ',  k_analog, ' +/- ', STD_analog

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
    call self % impProd % addScore(s1)
    call self % impAbs  % addScore(s2)

    ! Increase collision count
    call self % collCount % addScore(p % w)

  end subroutine reportInColl

  !!
  !! Process history report
  !! ASSUMPTIONS:
  !! **** FATE CODES NEED TO BE SPECIFIED
  !!
  subroutine reportHist(self,pre,post,fate)
    class(keffClerk), intent(inout) :: self
    class(phaseCoord), intent(in)        :: pre
    class(particle), intent(in)          :: post
    integer(shortInt),intent(in)         :: fate



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
    real(defReal)                        :: nuFiss, absorb, collCount, dummy, k_cycle

    ! Obtain end of cycle weight and k value used to change fission site generation rate
    endWgt = end % popWeight()
    k_cycle = end % k_eff

    ! Calculate and score analog estimate of k-eff
    k_est =  endWgt / self % startWgt * k_cycle

    !self % k_csum  = self % k_csum + k_est
    !self % k2_csum = self % k2_csum + k_est * k_est

    call self % k_analog % addScore(k_est)

    ! Calculate and score implicit estimate of k_eff
    call self % collCount % getScore(collCount,dummy,1)
    call self % impProd % getScore(nuFiss,dummy,collCount)
    call self % impAbs  % getScore(absorb,dummy,collCount)

    k_est = nuFiss/absorb

    call self % k_imp % addScore(k_est)

    call self % impProd % reset()
    call self % impAbs % reset()
    call self % collCount % reset()

    self % cycleCount = self % cycleCount + 1

  end subroutine reportCycleEnd

end module keffClerk_class
