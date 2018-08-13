module macroClerk_class_module

  use numPrecision
  use tallyCodes
  use genericProcedures, only : fatalError
  use dictionary_class,  only : dictionary
  use grid_class,        only : grid
  use particle_class,             only : particle, phaseCoord
  use particleDungeon_class,      only : particleDungeon
  use keffClerk_inter,            only : keffClerk
  use tallyEstimator_class,       only : tallyEstimator
  use transportNuclearData_inter, only : transportNuclearData
  use outputFile_class,           only : outputFile

  use xsMacroSet_class,           only : xsMacroSet_ptr

  implicit none
  private

  character(*),parameter :: CLASS_NAME = 'macroClerk'

  !interface macroClerk
  !  module procedure new_macroClerk
  !end interface

  !!
  !! Simple clerk to tally macroscopic reaction rates
  !! For now supports only collision estimator
  !!
  type, public :: macroClerk
    private
    character(nameLen)   :: name
    integer(shortInt)    :: cycleCount = 0                ! Cycles counter
    integer(shortInt)    :: N = 0                         ! Number of bins
    real(defReal)        :: targetRelSD = 0.0

    type(grid)                                    :: map  ! Map
    type(tallyEstimator),dimension(:),allocatable :: bins ! Result estimates

  contains
    ! Deferred Interface Procedures
    procedure :: validReports
    procedure :: display
    procedure :: isConverged
    procedure :: init
    !procedure :: print

    ! Overwrite report procedures
    !procedure :: reportInColl
    !procedure :: reportCycleEnd
  end type macroClerk

contains

  !!
  !! Return codes for reports this clerk accepts
  !!
  function validReports(self) result(validCodes)
    class(macroClerk), intent(in)              :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    validCodes = [inColl_CODE, cycleEnd_CODE]

  end function validReports

  !!
  !! Display progress. Print total estimate and maximum relative STD s
  !!
  subroutine display(self)
    class(macroClerk), intent(in)     :: self
    real(defReal),dimension(self % N) :: res
    real(defReal),dimension(self % N) :: std
    real(defReal)                     :: resSum, maxSTD

    ! Obtain current estimate of total rate and maximum std
    call self % bins % getEstimate(res, std, self % cycleCount)
    resSum = sum(res)

    ! Find maximum relative std
    where(res > FP_REL_TOL)
      std = std / res
    elsewhere
      std = ZERO
    end where
    maxSTD = maxval(std)

    ! Print estimates to a console
    print '(A,ES12.5,A,ES12.5)', 'Total rate: ', resSum, ' Maximum Rel. STD: ', maxSTD

  end subroutine display

  !!
  !! Perform convergance check in the Clerk
  !!
  function isConverged(self) result(isIt)
    class(macroClerk), intent(in)     :: self
    logical(defBool)                  :: isIt
    real(defReal),dimension(self % N) :: res
    real(defReal),dimension(self % N) :: std
    real(defReal)                     :: maxSTD

    ! Obtain current maximum STD
    call self % bins % getEstimate(res, std, self % cycleCount)

    ! Find maximum relative std
    where(res > FP_REL_TOL)
      std = std / res
    elsewhere
      std = ZERO
    end where
    maxSTD = maxval(std)

    ! Compare against target
    isIt = (maxSTD < self % targetRelSD)

  end function isConverged

  !!
  !! Initialise macroClerk from dictionary
  !!
  subroutine init(self,dict,name)
    class(macroClerk),intent(inout) :: self
    class(dictionary), intent(in)   :: dict
    character(nameLen), intent(in)  :: name
    character(nameLen)              :: type
    character(nameLen)              :: grid
    real(defReal)                   :: mini, maxi
    integer(shortInt)               :: N
    character(100),parameter :: Here ='init (keffActiveClerk_class.f90)'

!    ! Assign name
!    self % name = name
!
!    call dict % getOrDefault(type,'trigger','no')
!
!    ! Read convergance target
!    if( charCmp(type,'yes')) then
!      call dict % get(self % targetRelSD,'RelSDtarget')
!
!    end if

    ! Read bining settings *** Make them more usefull
!    if( dict % isPresent('grid')) then
!      call dict % get(grid,'grid')
!      call dict % get(mini,'min')
!      call dict % get(maxi,'max')
!      call dict % get(N,'bins')
!      call dict % get(type,'spacing')
!    else
!
!    end if

  end subroutine init


end module macroClerk_class_module
