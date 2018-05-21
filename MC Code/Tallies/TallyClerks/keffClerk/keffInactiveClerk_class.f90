module keffInactiveClerk_class

  use numPrecision
  use tallyCodes
  use particle_class,             only : particle, phaseCoord
  use particleDungeon_class,      only : particleDungeon
  use keffClerk_inter,            only : keffClerk
  use tallyEstimator_class,       only : tallyScore, tallyCounter

  implicit none
  private

  type, public,extends(keffClerk) :: keffInactiveClerk
    private
    real(defReal)             :: k_est
    real(defReal)             :: startWgt

  contains
    ! Deferred Interface Procedures
    procedure :: validReports
    procedure :: display

    ! Overwrite report procedures
    procedure :: reportCycleStart
    procedure :: reportCycleEnd

    ! Local procedures
    procedure :: keff

  end type keffInactiveClerk

contains

  !!
  !! Return codes for reports this clerk accepts
  !!
  function validReports(self) result(validCodes)
    class(keffInactiveClerk), intent(in)       :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    validCodes = [cycleStart_CODE, cycleEnd_CODE]

  end function validReports

  !!
  !! Display progress current estimate of k-eff with STD
  !!
  subroutine display(self)
    class(keffInactiveClerk), intent(in) :: self
    real(defReal)                        :: k_analog, STD_analog

    ! Obtain k-eff estimate for this cycle
    k_analog = self % k_est

    ! Print estimates to a console
    print '(A,F8.5,A,F8.5)', 'k-eff (analog): ',  k_analog

  end subroutine display

  !!
  !! Process beginning of a cycle
  !!
  subroutine reportCycleStart(self,start)
    class(keffInactiveClerk), intent(inout) :: self
    class(particleDungeon), intent(in)      :: start

    self % startWgt  = start % popWeight()

  end subroutine reportCycleStart

  !!
  !! Process end of the cycle
  !!
  subroutine reportCycleEnd(self,end)
    class(keffInactiveClerk), intent(inout) :: self
    class(particleDungeon), intent(in)      :: end
    real(defReal)                           :: endWgt, k_cycle

    ! Obtain end of cycle weight and k value used to change fission site generation rate
    endWgt = end % popWeight()
    k_cycle = end % k_eff

    ! Calculate and score analog estimate of k-eff
    self % k_est =  endWgt / self % startWgt * k_cycle

  end subroutine reportCycleEnd

  !!
  !! Return current estimate of k-eff
  !!
  pure function keff(self) result(k)
    class(keffInactiveClerk), intent(in) :: self
    real(defReal)                        :: k
    real(defReal)                        :: STD_analog

    k = self % k_est

  end function keff

end module keffInactiveClerk_class
