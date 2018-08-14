module tallyActiveAdmin_class

  use numPrecision
  use tallyCodes
  use genericProcedures,       only : charCmp
  use dictionary_class,        only : dictionary
  use tallyAdminBase_class,    only : tallyAdminBase, &
                                      reportInColl_super => reportInColl, &
                                      reportOutColl_super => reportOutColl, &
                                      reportHist_super => reportHist, &
                                      reportCycleStart_super => reportCycleStart, &
                                      reportCycleEnd_super => reportCycleEnd, &
                                      isConverged_super => isConverged, &
                                      init_super  => init ,&
                                      print_super => print, &
                                      kill_super  => kill ,&
                                      display_super => display
  use particle_class,          only : particle, phaseCoord
  use particleDungeon_class,   only : particleDungeon
  use keffActiveClerk_class,   only : keffActiveClerk
  use outputFile_class,        only : outputFile

  implicit none
  private

  type, public,extends(tallyAdminBase) :: tallyActiveAdmin
    private
    type(keffActiveClerk) :: keff_estimator
    logical(defBool)      :: keff_convergence = .false.
  contains
    ! New Interface
    procedure :: keff

    ! Extend superclass procedures
    procedure :: reportInColl
    procedure :: reportOutColl
    procedure :: reportHist
    procedure :: reportCycleStart
    procedure :: reportCycleEnd
    procedure :: init
    procedure :: print
    procedure :: kill
    procedure :: display
    procedure :: isConverged

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
  !! Process outgoing collision report
  !!
  subroutine reportOutColl(self,pre,post,MT,muL)
    class(tallyActiveAdmin), intent(inout) :: self
    class(phaseCoord), intent(in)          :: pre
    class(particle), intent(in)            :: post
    integer(shortInt), intent(in)          :: MT
    real(defReal), intent(in)              :: muL

    ! Process report with internal Clerk
    call self % keff_estimator % reportOutColl(pre, post, MT, muL)

    ! Call superclass procedure on self
    call reportOutColl_super(self, pre, post, MT, muL)

  end subroutine reportOutColl

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

    ! Process report with internal Clerk
    call self % keff_estimator % reportCycleEnd(end)

    ! Call superclass procedure on self
    call reportCycleEnd_super(self, end)

  end subroutine reportCycleEnd

  !!
  !! Initialise active admin
  !!
  subroutine init(self,dict)
    class(tallyActiveAdmin), intent(inout) :: self
    class(dictionary),intent(in)           :: dict
    type(dictionary)                       :: embDict
    character(nameLen)                     :: entry
    real(defReal)                          :: temp

    call embDict % init(3)

    ! Read settings for embedded clerk and put them into embDict
    call dict % getOrDefault(entry,'trigger','no')

    if( charCmp(entry,'yes') ) then
      self % keff_convergence = .true.
      call dict % get(temp,'SDtarget')
      call embDict % store('trigger','yes')
      call embDict % store('SDtarget',temp)
    end if

    call embDict % store('type','keffActiveClerk')

    ! Initialise embedded clerk
    entry = 'A_ADMIN'
    call self % keff_estimator % init(embDict,entry)

    ! Load rest of the clerks
    call init_super(self,dict)

  end subroutine init

  !!
  !! Add all results to outputfile
  !!
  subroutine print(self,output)
    class(tallyActiveAdmin), intent(in)  :: self
    class(outputFile), intent(inout)     :: output

    call self % keff_estimator % print(output)
    call print_super(self,output)

  end subroutine print

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

    print *,'Active Cycle:'
    call self % keff_estimator % display()

    ! Call superclass procedure on self
    call display_super(self)

  end subroutine display

  !!
  !!
  !!
  function isConverged(self) result(isIt)
    class(tallyActiveAdmin), intent(in)    :: self
    logical(defBool)                       :: isIt

    if( self % keff_convergence) then
      isIt = self % keff_estimator % isConverged()

      if ( self % checkConvergence ) isIt = isIt .and. isConverged_super(self)

    else if ( self % checkConvergence ) then
       isIt = isConverged_super(self)

    else
      isIt = .false.

    end if

  end function isConverged
    
end module tallyActiveAdmin_class
