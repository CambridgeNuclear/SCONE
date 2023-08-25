module tallyTimeAdmin_class

  use numPrecision
  use tallyCodes
  use genericProcedures,       only : charCmp
  use dictionary_class,        only : dictionary
  use tallyAdminBase_class,    only : tallyAdminBase, &
                                      reportInColl_super => reportInColl, &
                                      reportHist_super => reportHist, &
                                      reportCycleStart_super => reportCycleStart, &
                                      reportCycleEnd_super => reportCycleEnd, &
                                      isConverged_super => isConverged, &
                                      init_super => init ,&
                                      print_super => print, &
                                      kill_super => kill ,&
                                      display_super => display
  use particle_class,          only : particle, phaseCoord
  use particleDungeon_class,   only : particleDungeon
  use timeClerk_class,         only : timeClerk
  use outputFile_class,        only : outputFile

  implicit none
  private

  type, public,extends(tallyAdminBase) :: tallyTimeAdmin
    private
    type(timeClerk) :: power_estimator
    logical(defBool):: power_convergence = .false.
  contains
    ! New Interface
    procedure :: power
    procedure :: stepLength
    procedure :: incrementStep

    ! Extend superclass procedures
    procedure :: reportInColl
    procedure :: reportHist
    procedure :: reportCycleStart
    procedure :: reportCycleEnd
    procedure :: init
    procedure :: print
    procedure :: kill
    procedure :: display
    procedure :: isConverged

  end type tallyTimeAdmin

contains
  !!
  !! Return estimate of power to be used in normalisation and results
  !!
  function power(self) result(p)
    class(tallyTimeAdmin), intent(in) :: self
    real(defReal)                     :: p

    p = self % power_estimator % power()

  end function power

  !!
  !! Return timestep length given an index
  !!
  function stepLength(self,idx) result(dt)
    class(tallyTimeAdmin), intent(in) :: self
    integer(shortInt), intent(in)     :: idx
    real(defReal)                     :: dt

    dt = self % power_estimator % stepLength(idx)

  end function stepLength

  !!
  !! Increment step count
  !!
  subroutine incrementStep(self)
    class(tallyTimeAdmin), intent(inout) :: self
    call self % power_estimator % incrementStep()
  end subroutine incrementStep

  !!
  !! Process incoming collision report
  !!
  subroutine reportInColl(self, p)
    class(tallyTimeAdmin), intent(inout) :: self
    class(particle), intent(in)            :: p

    ! Process report with internal Clerk
    call self % power_estimator % reportInColl(p)

    ! Call superclass procedure on self
    call reportInColl_super(self, p)

  end subroutine reportInColl

  !!
  !! Process history report
  !! ASSUMPTIONS:
  !!
  subroutine reportHist(self, p)
    class(tallyTimeAdmin), intent(inout) :: self
    class(particle), intent(in)          :: p

    ! Process report with internal Clerk
    call self % power_estimator % reportHist(p)

    ! Call superclass procedure on self
    call reportHist_super(self, p)

  end subroutine reportHist

  !!
  !! Process beginning of a cycle
  !!
  subroutine reportCycleStart(self, start)
    class(tallyTimeAdmin), intent(inout) :: self
    class(particleDungeon), intent(in)   :: start

    ! Process report with internal Clerk
    call self % power_estimator % reportCycleStart(start)

    ! Call superclass procedure on self
    call reportCycleStart_super(self, start)

  end subroutine reportCycleStart

  !!
  !! Process end of the cycle
  !!
  subroutine reportCycleEnd(self, end)
    class(tallyTimeAdmin), intent(inout) :: self
    class(particleDungeon), intent(in)   :: end

    ! Process report with internal Clerk
    call self % power_estimator % reportCycleEnd(end)

    ! Call superclass procedure on self
    call reportCycleEnd_super(self, end)

  end subroutine reportCycleEnd

  !!
  !! Add all results to outputfile
  !!
  subroutine print(self,output)
    class(tallyTimeAdmin), intent(in)    :: self
    class(outputFile), intent(inout)     :: output

    print *, 'HERE'
    call self % power_estimator % print(output)
    call print_super(self,output)

  end subroutine print

  !!
  !! Initialise time admin
  !!
  subroutine init(self,dict)
    class(tallyTimeAdmin), intent(inout) :: self
    class(dictionary),intent(in)         :: dict
    character(nameLen)                   :: entry

    ! Initialise embedded clerk
    entry = 'T_ADMIN'
    call self % power_estimator % init(dict,entry)

    ! Load rest of the clerks
    call init_super(self,dict)

  end subroutine init

  !!
  !! Deallocates all content
  !!
  subroutine kill(self)
    class(tallyTimeAdmin), intent(inout) :: self

    ! Kill internal Clerks

    ! Call superclass procedure on self
    call kill_super(self)

  end subroutine kill

  !!
  !! Display results at the end of a cycle
  !!
  subroutine display(self)
    class(tallyTimeAdmin), intent(in) :: self

    call self % power_estimator % display()

    ! Call superclass procedure on self
    call display_super(self)

  end subroutine display

  !!
  !!
  !!
  function isConverged(self) result(isIt)
    class(tallyTimeAdmin), intent(in)    :: self
    logical(defBool)                     :: isIt

    if( self % power_convergence) then
      isIt = self % power_estimator % isConverged()

      if ( self % checkConvergence ) isIt = isIt .and. isConverged_super(self)

    else if ( self % checkConvergence ) then
       isIt = isConverged_super(self)

    else
      isIt = .false.

    end if

  end function isConverged

end module tallyTimeAdmin_class
