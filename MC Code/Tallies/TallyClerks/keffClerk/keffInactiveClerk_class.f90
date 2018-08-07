module keffInactiveClerk_class

  use numPrecision
  use tallyCodes
  use genericProcedures,          only : fatalError, charCmp
  use dictionary_class,           only : dictionary
  use particle_class,             only : particle, phaseCoord
  use particleDungeon_class,      only : particleDungeon
  use keffClerk_inter,            only : keffClerk
  use tallyEstimator_class,       only : tallyScore, tallyCounter

  implicit none
  private

  character(*),parameter :: CLASS_NAME = 'keffInactiveClerk'

  interface keffInactiveClerk
    module procedure new_keffInactiveClerk
  end interface


  type, public,extends(keffClerk) :: keffInactiveClerk
    private
    real(defReal)             :: k_est
    real(defReal)             :: startWgt

  contains
    ! Deferred Interface Procedures
    procedure :: validReports
    procedure :: display
    procedure :: init

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
    real(defReal)                        :: k_analog

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
  !! Initialise keffActiveClerk from dictionary
  !! Checks if type agrees with class name. if not returns error
  !!
  subroutine init(self,dict)
    class(keffInactiveClerk),intent(inout) :: self
    class(dictionary), intent(in)          :: dict
    character(nameLen)                     :: type
    character(100),parameter :: Here ='init (keffInctiveClerk_class.f90)'

    ! Check that class description matches class name
    call dict % get(type,'type')
!    if( .not.charCmp(CLASS_NAME, type) ) then
!!      call fatalError(Here, 'Type : ' // type // ' is different form class name ' // CLASS_NAME )

    if( dict % isPresent('trigger')) then
      call fatalError(Here,'keffInactiveClerk cannot be a convergance trigger')

    end if



  end subroutine init

  !!
  !! Return current estimate of k-eff
  !!
  pure function keff(self) result(k)
    class(keffInactiveClerk), intent(in) :: self
    real(defReal)                        :: k

    k = self % k_est

  end function keff

  !!
  !! keffInactiveClerk constructor
  !!
  function new_keffInactiveClerk(dict) result(new)
    class(dictionary),intent(in) :: dict
    type(keffInactiveClerk)      :: new

    call new % init(dict)

  end function new_keffInactiveClerk
end module keffInactiveClerk_class
