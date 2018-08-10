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
  use outputFile_class,           only : outputFile

  use xsMacroSet_class,           only : xsMacroSet_ptr

  implicit none
  private

  character(*),parameter :: CLASS_NAME = 'timeClerk'

  interface timeClerk
    module procedure new_timeClerk
  end interface

  type, public,extends(tallyClerk) :: timeClerk
    private
    character(nameLen)                               :: name             !! Clerk name

    integer(shortInt)                                :: stepCount = 1    !! Current step

    type(tallyCounter)                               :: impFission       !! Implicit fission
    type(tallyScore)                                 :: anaLeak          !! Analog neutron leakage

    integer(shortInt)                                :: nSteps           !! Number of timesteps
    real(defReal), dimension(:), allocatable, public :: stepLength       !! Time step length
    real(defReal), dimension(:), allocatable         :: power_imp        !! Stores power at each time step
    real(defReal), dimension(:), allocatable         :: power_std        !! Stores power STD at each time step
    real(defReal)                                    :: power0           !! Initial power
    real(defReal)                                    :: normFactor       !! Normalisation factor
    integer(shortInt)                                :: eventCounter = 0 !! Count collision events

    real(defReal)                                    :: startWgt
    real(defReal)                                    :: targetSD = 0.0

  contains
    ! Deferred Interface Procedures
    procedure :: validReports
    procedure :: display
    procedure :: isConverged
    procedure :: init
    procedure :: print

    ! Overwrite report procedures
    procedure :: reportInColl
    procedure :: reportHist
    procedure :: reportCycleStart
    procedure :: reportCycleEnd
    procedure :: power
    procedure :: incrementStep
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

    ! Print estimates to a console
    print '(A,ES15.5,A,ES15.5)', 'Power (implicit): ', &
    self % power_imp(self % stepCount), ' +/- ', self % power_std(self % stepCount)

  end subroutine display

  !!
  !! Perform convergance check in the Clerk
  !!
  !! NEVER USE THIS! Not Implemented yet.
  !!
  function isConverged(self) result(isIt)
    class(timeClerk), intent(in) :: self
    logical(defBool)             :: isIt
!    real(defReal)                :: power, SD

   ! SD = self % power_std(self % stepCount)

   ! isIt = (SD < self % targetSD)
    isIt = .false.

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
          call xsData % getMatMacroXS(XSs, p, p % matIdx())

        class default
          call fatalError(Here,'Dynamic type of XS data attached to particle is not transportNuclearData')

      end select
    end associate

    self % eventCounter = self % eventCounter + 1

    totalXS  = XSs % totalXS()
    fissXS   = XSs % fissionXS()

    ! Calculate flux and fission scores
    flux = p % w / totalXS

    s1 = fissXS * flux

    ! Add score to counter
    call self % impFission % addEstimate(s1)

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
    self % eventCounter = 0

    ! Reset score counters
    call self % impFission % reset()

  end subroutine reportCycleStart

  !!
  !! Process end of the time step
  !!
  subroutine reportCycleEnd(self,end)
    class(timeClerk), intent(inout)      :: self
    class(particleDungeon), intent(in)   :: end
    real(defReal)                        :: power_est, power_STD
    real(defReal)                        :: fission_imp, STD_imp

    ! Calculate implicit estimate of power
    call self % impFission % getEstimate(fission_imp, STD_imp, self % eventCounter)
    fission_imp = fission_imp * self % eventCounter
    STD_imp = STD_imp * self % eventCounter
    power_est = joulesPerMeV * energyPerFission * fission_imp / self % stepLength(self % stepCount)
    power_STD = joulesPerMeV * energyPerFission * STD_imp / self % stepLength(self % stepCount)

    ! If the first step, must normalise estimate to the initial power
    ! Normalisation factor must be stored for subsequent steps
    ! Don't require joulesPerMeV and energyPerFission should cancel - used here for future proofing
    ! when energyPerFission may be tallied explicitly
    if (self % stepCount == 1) then
      self % normFactor = self % power0 / power_est
      power_est = self % power0
      power_STD = power_STD * self % normFactor
    else
      power_est = power_est * self % normFactor
      power_STD = power_STD * self % normFactor
    end if

    self % power_imp(self % stepCount) = power_est
    self % power_std(self % stepCount) = power_STD

  end subroutine reportCycleEnd

  !!
  !! Write contents of the time Clerk to output file
  !!
  subroutine print(self,outFile)
    class(timeClerk), intent(in)       :: self
    class(outputFile), intent(inout)   :: outFile
    integer(shortInt)                  :: N,i
    real(defReal)                      :: time
    character(nameLen)                 :: name

    call outFile % startBlock(self % name)

    ! Calculate grid size
    N = size(self % stepLength)

    ! Print time grid. Assume start at time = 0
    name = 'timeGrid_ends'
    time = ZERO
    call outFile % startArray(name,[N])

    do i=1,N
      time = time + self % stepLength(i)
      call outFile % addValue(time)
    end do

    call outFile % endArray()

    ! Print power results
    name = 'power'
    call outFile % startArray(name,[N])

    call outFile % addResult(self % power_imp, self % power_std)

    call outFile % endArray()

    call outFile % endBlock()

  end subroutine print


  !!
  !! Initialise timeClerk from dictionary
  !! Checks if type agrees with class name. if not returns error
  !!
  subroutine init(self,dict, name)
    class(timeClerk),intent(inout)       :: self
    class(dictionary), intent(in)        :: dict
    character(nameLen),intent(in)        :: name
    character(nameLen)                   :: type
    character(100),parameter :: Here ='init (timeClerk_class.f90)'

    ! Assign name
    self % name = name

    call dict % getOrDefault(type,'trigger','no')

    ! Read convergence target
    if( charCmp(type,'yes')) then
      call dict % get(self % targetSD,'SDtarget')

    end if

    ! Read initial power
    call dict % get(self % power0, 'power')

    ! Read number of time steps
    call dict % get(self % nSteps, 'nsteps')
    allocate(self % power_imp(self % nSteps))
    allocate(self % power_std(self % nSteps))

    ! Read array providing time step lengths
    call dict % get(self % stepLength, 'dt')

  end subroutine init

  !!
  !! Return estimate of power
  !!
  function power(self) result(p)
    class(timeClerk), intent(in) :: self
    real(defReal)                :: p

    ! On the first step, use the specified initial power
    if (self % stepCount == 1) then
      p = self % power0
    else
      p = self % power_imp(self % stepCount)
    end if

  end function power

  !!
  !! Increment step count
  !!
  subroutine incrementStep(self)
    class(timeClerk), intent(inout) :: self
    self % stepCount = self % stepCount + 1
  end subroutine incrementStep

  !!
  !! timeClerk constructor function
  !!
  function new_timeClerk(dict,name) result(new)
    class(dictionary), intent(in)  :: dict
    character(nameLen), intent(in) :: name
    type(timeClerk)                :: new

    call new % init(dict,name)

  end function new_timeClerk

end module timeClerk_class
