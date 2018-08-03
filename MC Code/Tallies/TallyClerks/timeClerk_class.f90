module timeClerk_class

  use numPrecision
  use universalVariables,         only : energyPerFiss, joulesPerMeV
  use tallyCodes
  use genericProcedures,          only : fatalError, charCmp
  use dictionary_class,           only : dictionary
  use particle_class,             only : particle, phaseCoord
  use particleDungeon_class,      only : particleDungeon
  use tallyClerk_inter,           only : tallyClerk
  use tallyEstimator_class,       only : tallyScore, tallyCounter
  use transportNuclearData_inter, only : transportNuclearData

  use xsMacroSet_class,           only : xsMacroSet_ptr

  implicit none
  private

  character(*),parameter :: CLASS_NAME = 'timeClerk'

  type, public,extends(tallyClerk) :: timeClerk
    private

    integer(shortInt)                    :: cycleCount = 0

    type(tallyScore)                     :: impFissions ! Implicit fission
    type(tallyScore)                     :: anaLeak     ! Analog neutron leakage

    type(tallyCounter)                   :: power_imp

    real(defReal)                        :: startWgt
    real(defReal)                        :: targetSD = 0.0

    real(defReal)                        :: stepLength = 1.0 ! Time step length

  contains
    ! Deferred Interface Procedures
    procedure :: validReports
    procedure :: display
    procedure :: isConverged
    procedure :: init


    ! Overwrite report procedures
    procedure :: reportInColl
    procedure :: reportHist
    procedure :: reportCycleStart
    procedure :: reportCycleEnd
    procedure :: power
  end type timeClerk

contains

  !!
  !! Return codes for reports this clerk accepts
  !!
  function validReports(self) result(validCodes)
    class(timeClerk), intent(in)         :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    validCodes = [inColl_CODE, cycleStart_CODE, cycleEnd_CODE, hist_CODE]

  end function validReports

  !!
  !! Display progress current estimate of k-eff with STD
  !!
  subroutine display(self)
    class(timeClerk), intent(in) :: self
    real(defReal)                :: power_imp, STD_imp, STD_analog

    ! Obtain current implicit estimate of power
    call self % power_imp % getEstimate(power_imp, STD_imp, self % stepCount)

    ! Print estimates to a console
    print '(A,F8.5,A,F8.5)', 'Power (implicit): ', power_imp, ' +/- ', STD_imp

  end subroutine display

  !!
  !! Perform convergance check in the Clerk
  !!
  function isConverged(self) result(isIt)
    class(timeClerk), intent(in) :: self
    logical(defBool)             :: isIt
    real(defReal)                :: power, SD

    call self % power_imp % getEstimate(power,SD,self % cycleCount)

    isIt = (SD < self % targetSD)

  end function isConverged


  !!
  !! Process collision report
  !!
  subroutine reportInColl(self,p)
    class(timeClerk), intent(inout) :: self
    class(particle), intent(in)           :: p
    type(xsMacroSet_ptr)                  :: XSs
    real(defReal)                         :: totalXS, fissXS, flux
    real(defReal)                         :: s1
    character(100), parameter  :: Here = 'reportInColl (timeClerk_class.f90)'

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
    fissXS   = XSs % fissionXS()

    ! Calculate flux and energy scores
    flux = p % w / totalXS

    s1 = fissXS * flux

    ! Add score to counter
    call self % impFission % add(s1)

  end subroutine reportInColl

  !!
  !! Process history report
  !!
  subroutine reportHist(self,pre,post,fate)
    class(timeClerk), intent(inout):: self
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
  !! Process beginning of a time step
  !!
  subroutine reportCycleStart(self,start)
    class(timeClerk), intent(inout)      :: self
    class(particleDungeon), intent(in)   :: start

    self % startWgt = start % popWeight()

  end subroutine reportCycleStart

  !!
  !! Process end of the time step
  !!
  subroutine reportCycleEnd(self,end)
    class(timeClerk), intent(inout)      :: self
    class(particleDungeon), intent(in)   :: end
    real(defReal)                        :: power_est
    real(defReal)                        :: fissions

    ! Obtain end of cycle weight and k value used to change fission site generation rate
    !endWgt = end % popWeight()
    !power_step = end % power

    ! Calculate and score analog estimate of power
    !power_est =  endWgt / self % startWgt * power_step
    !call self % power_analog % addEstimate(power_est)

    ! Calculate and score implicit estimate of power
    fissions  = self % impFission % get()

    power_est = joulesPerMeV * energyPerFission * fissions / self % stepLength

    call self % power_imp % addEstimate(power_est)

    ! Reset score counters
    call self % impFission % reset()

    ! Increas counter of steps
    self % stepCount = self % stepCount + 1

  end subroutine reportCycleEnd

  !!
  !! Initialise timeClerk from dictionary
  !! Checks if type agrees with class name. if not returns error
  !!
  subroutine init(self,dict)
    class(timeClerk),intent(inout) :: self
    class(dictionary), intent(in)        :: dict
    character(nameLen)                   :: type
    character(100),parameter :: Here ='init (timeClerk_class.f90)'

    call dict % get(type,'type')

    ! Check that class description matches class name
    if ( .not.charCmp(CLASS_NAME,type) ) then
      call fatalError(Here, 'Type : ' // type // ' is different form class name ' // CLASS_NAME )

    end if

    call dict % getOrDefault(type,'trigger','no')

    ! Read convergence target
    if( charCmp(type,'yes')) then
      call dict % get(self % targetSD,'SDtarget')

    end if

  end subroutine init


  !!
  !! Return estimate of power
  !!
  function power(self) result(p)
    class(keffActiveClerk), intent(in) :: self
    real(defReal)                      :: p

    call self % power_imp % getEstimate(p, self % stepCount)

  end function power

  !!
  !! timeClerk constructor function
  !!
  function new_timeClerk(dict) result(new)
    class(dictionary), intent(in) :: dict
    type(timeClerk)               :: new

    call new % init(dict)

  end function new_timeClerk

end module timeClerk_class
