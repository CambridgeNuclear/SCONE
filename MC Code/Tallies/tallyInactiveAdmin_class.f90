module tallyInactiveAdmin_class

  use numPrecision
  use tallyCodes
  use tallyAdminBase_class,    only : tallyAdminBase, &
                                      reportCycleStart_super => reportCycleStart, &
                                      reportCycleEnd_super => reportCycleEnd, &
                                      kill_super => kill
  use particle_class,          only : particle, phaseCoord
  use particleDungeon_class,   only : particleDungeon
  use keffInactiveClerk_class, only : keffInactiveClerk

  implicit none
  private

  type, public,extends(tallyAdminBase) :: tallyInactiveAdmin
    private
    type(keffInactiveClerk) :: keff_estimator
  contains
    ! New interface
    procedure :: keff

    ! Extend superclass procedures
    procedure :: reportCycleStart
    procedure :: reportCycleEnd
    procedure :: kill

    ! Overwrite superclass procedures
    procedure :: display


  end type tallyInactiveAdmin

contains

  !!
  !! Return estimate of k-eff to be used in a calculation for normalisation
  !!
  function keff(self) result(k)
    class(tallyInactiveAdmin), intent(in) :: self
    real(defReal)                         :: k

    k = self % keff_estimator % keff()

  end function keff

  !!
  !! Process beginning of a cycle
  !!
  subroutine reportCycleStart(self,start)
    class(tallyInactiveAdmin), intent(inout) :: self
    class(particleDungeon), intent(in)       :: start
    integer(shortInt)                        :: i, idx

    ! Process report with internal Clerk
    call self % keff_estimator % reportCycleStart(start)

    ! Call superclass procedure on self
    call reportCycleStart_super(self,start)

  end subroutine reportCycleStart

  !!
  !! Process end of the cycle
  !!
  subroutine reportCycleEnd(self,end)
    class(tallyInactiveAdmin), intent(inout) :: self
    class(particleDungeon), intent(in)       :: end
    integer(shortInt)                        :: i, idx

    ! Process report with internal Clerk
    call self % keff_estimator % reportCycleEnd(end)

    ! Call superclass procedure on self
    call reportCycleEnd_super(self,end)

  end subroutine reportCycleEnd

  !!
  !! Deallocates all content
  !!
  subroutine kill(self)
    class(tallyInactiveAdmin), intent(inout) :: self

    ! Kill internal Clerks

    ! Call superclass procedure on self
    call kill_super(self)

  end subroutine kill

  !!
  !! Display results at the end of a cycle
  !!
  subroutine display(self)
    class(tallyInactiveAdmin), intent(in) :: self

    print *,'Inactive Cycle:'
    call self % keff_estimator % display()

  end subroutine display
    
end module tallyInactiveAdmin_class
