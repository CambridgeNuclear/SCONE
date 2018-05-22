module tallyInactiveAdmin_class

  use numPrecision
  use tallyCodes
  use tallyAdminBase_class,    only : tallyAdminBase, &
                                      reportCycleStartBase => reportCycleStart, &
                                      reportCycleEndBase => reportCycleEnd, &
                                      killBase => kill
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

    ! Extend Base class procedures
    procedure :: reportCycleStart
    procedure :: reportCycleEnd
    procedure :: kill

    ! Overwrite Base class procedures
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

    ! Call base class procedure on self
    call reportCycleStartBase(self,start)

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

    ! Call base class procedure on self
    call reportCycleEndBase(self,end)

  end subroutine reportCycleEnd

  !!
  !! Deallocates all content
  !!
  subroutine kill(self)
    class(tallyInactiveAdmin), intent(inout) :: self

    ! Kill internal Clerks

    ! Call base class procedure on self
    call killBase(self)

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
