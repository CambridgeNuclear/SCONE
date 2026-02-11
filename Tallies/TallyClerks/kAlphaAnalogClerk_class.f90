module kAlphaAnalogClerk_class

  use numPrecision
  use tallyCodes
  use dictionary_class,      only : dictionary
  use genericProcedures,     only : fatalError
  use particle_class,        only : particle
  use particleDungeon_class, only : particleDungeon
  use outputFile_class,      only : outputFile

  use scoreMemory_class,     only : scoreMemory
  use tallyResult_class,     only : tallyResult, tallyResultEmpty
  use tallyClerk_inter,      only : tallyClerk, kill_super => kill

  implicit none
  private

  !!
  !! Simplest possible analog k estimator that determines
  !! criticality by comparing population weight before and after transport cycle.
  !! Uses this to update an estimate of the alpha eigenvalue.
  !!
  !! The update depends on whether the system is believed to be sub- or
  !! super-critical. Will give dubious results if the calculation is initialised
  !! in the wrong k-regime.
  !!
  !! Note that this estimator suffers from large inter-cycle correlations and
  !! tends to have an unreliable estimate of uncertainty.
  !!
  !! Private Members:
  !!   startPopWgt -> Accumulated total particle Starting weight from cycle
  !!   endPopWgt   -> Accumulated total particle End weight from previous cycle
  !!   alpha       -> Running estimate of alpha
  !!   k           -> Running estimate of k
  !!
  !! SAMPLE DICTIOANRY INPUT:
  !!
  !! myClerk {
  !!   type kAlphaAnalogClerk;
  !! }
  !!
  type, public,extends(tallyClerk) :: kAlphaAnalogClerk
    private
    real(defReal) :: startPopWgt = ZERO
    real(defReal) :: endPopWgt   = ZERO
    real(defReal) :: k           = ONE
    real(defReal) :: alpha       = ZERO
  contains
    ! Procedures used during build
    procedure :: init
    procedure :: kill
    procedure :: validReports
    procedure :: getSize

    ! File reports and check status -> run-time procedures
    procedure :: reportCycleStart
    procedure :: reportCycleEnd

    ! Output procedures
    procedure  :: display
    procedure  :: print
    procedure  :: getResult
  end type kAlphaAnalogClerk

  !!
  !! k/alpha result class
  !!
  !! Public Members:
  !!   k     -> Current k estimate
  !!   alpha -> Current alpha estimate
  !!
  type,public, extends(tallyResult) :: kAlphaResult
    real(defReal) :: k = ONE
    real(defReal) :: alpha = ZERO
  contains
    procedure :: initRes
  end type kAlphaResult


contains

  !!
  !! Initialise from dictionary and name
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine init(self, dict, name)
    class(kAlphaAnalogClerk), intent(inout) :: self
    class(dictionary), intent(in)        :: dict
    character(nameLen), intent(in)       :: name
    character(100),parameter :: Here = 'init (kAlphaAnalogClerk.f90)'

    ! Needs no settings, just load name
    call self % setName(name)

    ! Ensure correct initialisation to default values
    self % startPopWgt = ZERO
    self % endPopWgt   = ZERO
    self % k = ONE

    ! Initial value for alpha
    call dict % getOrDefault(self % alpha, 'alpha_0', 0.1_defReal)

  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(kAlphaAnalogClerk), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    self % startPopWgt = ZERO
    self % endPopWgt = ZERO
    self % alpha = ZERO
    self % k = ONE

  end subroutine kill

  !!
  !! Returns array of codes that represent diffrent reports
  !!
  !! See tallyClerk_inter for details
  !!
  function validReports(self) result(validCodes)
    class(kAlphaAnalogClerk),intent(in)        :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    validCodes = [ cycleStart_CODE, cycleEnd_CODE ]

  end function validReports

  !!
  !! Return memory size of the clerk
  !!
  !! See tallyClerk_inter for details
  !!
  elemental function getSize(self) result(S)
    class(kAlphaAnalogClerk), intent(in) :: self
    integer(shortInt)                    :: S

    S = 2

  end function getSize

  !!
  !! Process beginning of a cycle
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportCycleStart(self, start, mem)
    class(kAlphaAnalogClerk), intent(inout) :: self
    class(particleDungeon), intent(in)      :: start
    type(scoreMemory), intent(inout)        :: mem

    ! Update start population weight
    self % startPopWgt = self % startPopWgt + start % popWeight()

  end subroutine reportCycleStart

  !!
  !! Process end of the cycle
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportCycleEnd(self, end, mem)
    class(kAlphaAnalogClerk), intent(inout) :: self
    class(particleDungeon), intent(in)      :: end
    type(scoreMemory), intent(inout)        :: mem
    real(defReal)                           :: k_norm

    ! Update end population weight
    self % endPopWgt = self % endPopWgt + end % popWeight()

    ! Close batch
    if( mem % lastCycle() ) then
      k_norm = end % k_eff
      
      ! Calculate and score analog estimate of k
      self % k = self % endPopWgt / self % startPopWgt * k_norm
      
      ! Update the estimate of alpha
      if (self % alpha >= ZERO) then
        self % alpha = self % alpha * self % k
      else
        self % alpha = self % alpha / self % k
      end if
      call mem % accumulate(self % k, self % getMemAddress() )
      call mem % accumulate(self % alpha, self % getMemAddress() + 1)

      self % startPopWgt = ZERO
      self % endPopWgt   = ZERO
    end if

  end subroutine reportCycleEnd

  !!
  !! Display convergance progress on the console
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine display(self, mem)
    class(kAlphaAnalogClerk), intent(in) :: self
    type(scoreMemory), intent(in)        :: mem
    real(defReal)                        :: k, STDk, alpha, STDa

    call mem % getResult(k, STDk, self % getMemAddress())
    call mem % getResult(alpha, STDa, self % getMemAddress() + 1)

    ! Print estimates to a console
    print '(A,F8.5,A,F8.5)', 'k (analog): ',  k, ' +/- ', STDk
    print '(A,F8.5)', 'k current: ',  self % k
    print '(A,F8.4,A,F8.4)', 'Alpha accumulated: ',  alpha/10**6, 'us^-1 +/- ', STDa/10**6
    print '(A,F8.4,A)', 'Alpha current: ',  self % alpha/10**6, 'us^-1'

  end subroutine display

  !!
  !! Write contents of the clerk in the slot to output file
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine print(self, outFile, mem)
    class(kAlphaAnalogClerk), intent(in) :: self
    class(outputFile), intent(inout)     :: outFile
    type(scoreMemory), intent(in)        :: mem
    real(defReal)                        :: k, STDk, alpha, STDa
    character(nameLen)                   :: name

    ! Get result value
    call mem % getResult(k, STDk, self % getMemAddress())
    call mem % getResult(alpha, STDa, self % getMemAddress() + 1)

    ! Print to output file
    call outFile % startBlock(self % getName() )
    name = 'k'
    call outFile % printResult(k, STDk, name)
    name = 'alpha'
    call outFile % printResult(alpha, STDa, name)
    call outFile % endBlock()

  end subroutine print

  !!
  !! Return result for interaction with Physics Package
  !! from the clerk in the slot
  !!
  !! See tallyClerk_inter for details
  !!
  pure subroutine getResult(self, res, mem)
    class(kAlphaAnalogClerk), intent(in)            :: self
    class(tallyResult), allocatable, intent(inout)  :: res
    type(scoreMemory), intent(in)                   :: mem
    type(kAlphaResult)                              :: myKAlpha

    ! Get result value
    call myKAlpha % initRes(self % k, self % alpha)
    allocate(res, source = myKAlpha)

  end subroutine getResult

  !!
  !! Initialise the kAlphaResult
  !!
  pure subroutine initRes(self, k, alpha)
    class(kAlphaResult), intent(inout) :: self
    real(defReal), intent(in)         :: k
    real(defReal), intent(in)         :: alpha

    self % k = k
    self % alpha = alpha

  end subroutine initRes

end module kAlphaAnalogClerk_class
