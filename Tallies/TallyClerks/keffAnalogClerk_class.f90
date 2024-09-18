module keffAnalogClerk_class

  use numPrecision
  use tallyCodes
  use universalVariables,    only : MAX_COL
  use dictionary_class,      only : dictionary
  use genericProcedures,     only : fatalError
  use display_func,          only : statusMsg
  use particle_class,        only : particle
  use particleDungeon_class, only : particleDungeon
  use outputFile_class,      only : outputFile

  use scoreMemory_class,     only : scoreMemory
  use tallyResult_class,     only : tallyResult, tallyResultEmpty
  use tallyClerk_inter,      only : tallyClerk, kill_super => kill

  implicit none
  private

  !! Locations of diffrent bins wrt memory address of the clerk
  integer(shortInt), parameter :: MEM_SIZE = 3
  integer(longInt), parameter  :: START_POP = 0 ,&  ! Population tally at the start of the cycle
                                  END_POP   = 1 ,&  ! Population tally at the end of the cycle
                                  K_EST     = 2     ! k-eff estimate

  !!
  !! Simplest possible analog k-eff estimator that determines
  !! criticality by comparing population weight before and after transport cycle
  !! Assumes that constant to normalise fission site production is uniform in geometry
  !!
  !! Private Members:
  !!   startPopWgt -> Accululated total particle Statring weight over ALL previous cycles
  !!   endPopWgt   -> Accululated total particle End weight over ALL previous cycles
  !!
  !! SAMPLE DICTIOANRY INPUT:
  !!
  !! myClerk {
  !!   type keffAnalogClerk;
  !! }
  !!
  type, public,extends(tallyClerk) :: keffAnalogClerk
  contains
    ! Procedures used during build
    procedure :: init
    procedure :: kill
    procedure :: validReports
    procedure :: getSize

    ! File reports and check status -> run-time procedures
    procedure :: reportCycleStart
    procedure :: reportCycleEnd
    procedure :: closeCycle

    ! Output procedures
    procedure  :: display
    procedure  :: print
    procedure  :: getResult

  end type keffAnalogClerk

  !!
  !! k-eff result class
  !!
  !! Public Members:
  !!   keff -> Result, keff(1) is criticality, keff(2) is STD
  !!
  type, public, extends(tallyResult) :: keffResult
    real(defReal), dimension(2) :: keff = [ONE, ZERO]
  end type keffResult


contains

  !!
  !! Initialise from dictionary and name
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine init(self, dict, name)
    class(keffAnalogClerk), intent(inout) :: self
    class(dictionary), intent(in)        :: dict
    character(nameLen), intent(in)       :: name
    character(100),parameter :: Here = 'init (keffAnalogClerk.f90)'

    ! Needs no settings, just load name
    call self % setName(name)

  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(keffAnalogClerk), intent(inout) :: self

    ! Superclass
    call kill_super(self)

  end subroutine kill

  !!
  !! Returns array of codes that represent diffrent reports
  !!
  !! See tallyClerk_inter for details
  !!
  function validReports(self) result(validCodes)
    class(keffAnalogClerk),intent(in)           :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    validCodes = [cycleStart_CODE, cycleEnd_CODE, closeCycle_CODE]

  end function validReports

  !!
  !! Return memory size of the clerk
  !!
  !! See tallyClerk_inter for details
  !!
  elemental function getSize(self) result(S)
    class(keffAnalogClerk), intent(in) :: self
    integer(shortInt)                  :: S

    S = MEM_SIZE

  end function getSize

  !!
  !! Process beginning of a cycle
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportCycleStart(self, start, mem)
    class(keffAnalogClerk), intent(inout) :: self
    class(particleDungeon), intent(in)    :: start
    type(scoreMemory), intent(inout)      :: mem

    ! Add score to counter
    call mem % score(start % popWeight(), self % getMemAddress() + START_POP)

  end subroutine reportCycleStart

  !!
  !! Process end of the cycle
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportCycleEnd(self, end, mem)
    class(keffAnalogClerk), intent(inout) :: self
    class(particleDungeon), intent(in)    :: end
    type(scoreMemory), intent(inout)      :: mem

    ! Add score to counter
    call mem % score(end % popWeight(), self % getMemAddress() + END_POP)

  end subroutine reportCycleEnd

  !!
  !! Close the cycle
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine closeCycle(self, end, mem)
    class(keffAnalogClerk), intent(inout) :: self
    class(particleDungeon), intent(in)    :: end
    type(scoreMemory), intent(inout)      :: mem
    real(defReal)                         :: k_norm, k_eff, startPopWgt, endPopWgt

    ! Close batch
    if (mem % lastCycle()) then
      k_norm = end % k_eff

      startPopWgt = mem % getScore(self % getMemAddress() + START_POP)
      endPopWgt   = mem % getScore(self % getMemAddress() + END_POP)

      ! Calculate and score analog estimate of k-eff
      k_eff =  endPopWgt / startPopWgt * k_norm
      call mem % accumulate(k_eff, self % getMemAddress() + K_EST)

    end if

  end subroutine closeCycle

  !!
  !! Display convergance progress on the console
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine display(self, mem)
    class(keffAnalogClerk), intent(in)  :: self
    type(scoreMemory), intent(in)       :: mem
    real(defReal)                       :: k, STD
    character(MAX_COL)                  :: buffer

    call mem % getResult(k, STD, self % getMemAddress() + K_EST)

    ! Print estimates to a console
    write(buffer, '(A,F8.5,A,F8.5)') 'k-eff (analog): ',  k, ' +/- ', STD
    call statusMsg(buffer)

  end subroutine display

  !!
  !! Write contents of the clerk in the slot to output file
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine print(self, outFile, mem)
    class(keffAnalogClerk), intent(in) :: self
    class(outputFile), intent(inout)   :: outFile
    type(scoreMemory), intent(in)      :: mem
    real(defReal)                      :: k, STD
    character(nameLen)                 :: name

    ! Get result value
    call mem % getResult(k, STD, self % getMemAddress() + K_EST)

    ! Print to output file
    call outFile % startBlock(self % getName())
    name = 'k_analog'
    call outFile % printResult(k, STD, name)
    call outFile % endBlock()

  end subroutine print

  !!
  !! Return result for interaction with Physics Package
  !! from the clerk in the slot
  !!
  !! See tallyClerk_inter for details
  !!
  pure subroutine getResult(self, res, mem)
    class(keffAnalogClerk), intent(in)              :: self
    class(tallyResult), allocatable, intent(inout)  :: res
    type(scoreMemory), intent(in)                   :: mem
    real(defReal)                                   :: k, STD

    ! Get result value
    call mem % getResult(k, STD, self % getMemAddress() + K_EST)

    allocate(res, source = keffResult([k, STD]))

  end subroutine getResult


end module keffAnalogClerk_class
