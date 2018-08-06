module timeClerk_class

  use numPrecision
  use universalVariables
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

  interface timeClerk
    module procedure new_timeClerk
  end interface

  type, public,extends(tallyClerk) :: timeClerk
    private

    integer(shortInt)                                :: stepCount = 1 ! Current step

    type(tallyScore)                                 :: impFission    ! Implicit fission
    type(tallyScore)                                 :: anaLeak       ! Analog neutron leakage

    integer(shortInt)                                :: nSteps        ! Number of timesteps
    real(defReal), dimension(:), allocatable, public :: stepLength    ! Time step length
    type(tallyCounter), dimension(:), allocatable    :: power_imp     ! Stores power at each time step
    real(defReal)                                    :: power0        ! Initial power
    real(defReal)                                    :: normFactor    ! Normalisation factor

    real(defReal)                                    :: startWgt
    real(defReal)                                    :: targetSD = 0.0

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
  !! Display progress current estimate of power with STD
  !!
  subroutine display(self)
    class(timeClerk), intent(in) :: self
    real(defReal)                :: power_imp, STD_imp, STD_analog

    ! Obtain current implicit estimate of power
    call self % power_imp(self % stepCount) % getEstimate(power_imp, STD_imp, 1)

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

    call self % power_imp(self % stepCount) % getEstimate(power,SD,1)

    isIt = (SD < self % targetSD)

  end function isConverged


  !!
  !! Process collision report
  !!
  subroutine reportInColl(self,p)
    class(timeClerk), intent(inout) :: self
    class(particle), intent(in)     :: p
    type(xsMacroSet_ptr)            :: XSs
    real(defReal)                   :: totalXS, fissXS, flux
    real(defReal)                   :: s1
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

    ! Calculate flux and fission scores
    flux = p % w / totalXS

    s1 = fissXS * flux

    ! Add score to counter
    call self % impFission % add(s1)

  end subroutine reportInColl

  !!
  !! Process history report
  !! Probably needn't be here!
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

    ! Calculate implicit estimate of power
    fissions  = self % impFission % get()

    power_est = joulesPerMeV * energyPerFission * fissions / self % stepLength(self % stepCount)

    ! If the first step, must normalise estimate to the initial power
    ! Normalisation factor must be stored for subsequent steps
    ! Don't require joulesPerMeV and energyPerFission should cancel - used here for future proofing
    ! when energyPerFission may be tallied explicitly
    if (self % stepCount == 1) then
      self % normFactor = self % power0 / power_est
      power_est = self % power0
    else
      power_est = power_est * self % normFactor
    end if

    call self % power_imp(self % stepCount) % addEstimate(power_est)

    ! Reset score counters
    call self % impFission % reset()

    ! Increment step count
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

!    call dict % get(type,'type')

!    ! Check that class description matches class name
!    if ( .not.charCmp(CLASS_NAME,type) ) then
!      call fatalError(Here, 'Type : ' // type // ' is different form class name ' // CLASS_NAME )
!
!    end if

    call dict % getOrDefault(type,'trigger','no')

    ! Read convergence target
    if( charCmp(type,'yes')) then
      call dict % get(self % targetSD,'SDtarget')

    end if

    ! Read initial power
    call dict % get(self % power0, 'power')

    ! Read number of time steps
    call dict % get(self % nSteps, 'numSteps')

    ! Read array providing time step lengths
    call dict % get(self % stepLength, 'steps')

  end subroutine init

  !!
  !! Return estimate of power
  !!
  function power(self) result(p)
    class(timeClerk), intent(in) :: self
    real(defReal)                :: p

    call self % power_imp(self % stepCount) % getEstimate(p, 1)

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
