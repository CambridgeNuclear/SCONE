module tallyActiveAdmin_class

  use numPrecision
  use tallyCodes
  use tallyAdminBase_class,    only : tallyAdminBase, &
                                      reportInColl_super => reportInColl, &
                                      reportHist_super => reportHist, &
                                      reportCycleStart_super => reportCycleStart, &
                                      reportCycleEnd_super => reportCycleEnd, &
                                      kill_super => kill
  use particle_class,          only : particle, phaseCoord
  use particleDungeon_class,   only : particleDungeon
  use keffActiveClerk_class,   only : keffActiveClerk

  implicit none
  private

  type, public,extends(tallyAdminBase) :: tallyActiveAdmin
    private
    type(keffActiveClerk) :: keff_estimator
  contains
    ! New Interface
    procedure :: keff

    ! Extend superclass procedures
    procedure :: reportInColl
    procedure :: reportHist
    procedure :: reportCycleStart
    procedure :: reportCycleEnd
    procedure :: kill

    ! Overwrite superclass procedures
    procedure :: display

  end type tallyActiveAdmin

contains
  !!
  !! Return estimate of k-eff to be used in a calculation for normalisation
  !!
  function keff(self) result(k)
    class(tallyActiveAdmin), intent(in) :: self
    real(defReal)                         :: k

    k = self % keff_estimator % keff()

  end function keff

  !!
  !! Process incoming collision report
  !!
  subroutine reportInColl(self, p)
    class(tallyActiveAdmin), intent(inout) :: self
    class(particle), intent(in)            :: p

    ! Process report with internal Clerk
    call self % keff_estimator % reportInColl(p)

    ! Call superclass procedure on self
    call reportInColl_super(self, p)

  end subroutine reportInColl

  !!
  !! Process history report
  !! ASSUMPTIONS:
  !!
  subroutine reportHist(self, pre, post, fate)
    class(tallyActiveAdmin), intent(inout) :: self
    class(phaseCoord), intent(in)          :: pre
    class(particle), intent(in)            :: post
    integer(shortInt),intent(in)           :: fate

    ! Process report with internal Clerk
    call self % keff_estimator % reportHist(pre, post, fate)

    ! Call superclass procedure on self
    call reportHist_super(self, pre, post, fate)

  end subroutine reportHist

  !!
  !! Process beginning of a cycle
  !!
  subroutine reportCycleStart(self, start)
    class(tallyActiveAdmin), intent(inout) :: self
    class(particleDungeon), intent(in)     :: start
    integer(shortInt)                      :: i, idx

    ! Process report with internal Clerk
    call self % keff_estimator % reportCycleStart(start)

    ! Call superclass procedure on self
    call reportCycleStart_super(self, start)

  end subroutine reportCycleStart

  !!
  !! Process end of the cycle
  !!
  subroutine reportCycleEnd(self, end)
    class(tallyActiveAdmin), intent(inout) :: self
    class(particleDungeon), intent(in)     :: end
    integer(shortInt)                      :: i, idx

    ! Process report with internal Clerk
    call self % keff_estimator % reportCycleEnd(end)

    ! Call superclass procedure on self
    call reportCycleEnd_super(self, end)

  end subroutine reportCycleEnd

  !!
  !! Deallocates all content
  !!
  subroutine kill(self)
    class(tallyActiveAdmin), intent(inout) :: self

    ! Kill internal Clerks

    ! Call superclass procedure on self
    call kill_super(self)

  end subroutine kill

  !!
  !! Display results at the end of a cycle
  !!
  subroutine display(self)
    class(tallyActiveAdmin), intent(in) :: self

    print *,'Inactive Cycle:'
    call self % keff_estimator % display()

  end subroutine display

    
end module tallyActiveAdmin_class
