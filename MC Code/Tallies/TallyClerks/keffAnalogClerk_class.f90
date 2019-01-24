module keffAnalogClerk_class

  use numPrecision
  use tallyCodes
  use dictionary_class,      only : dictionary
  use genericProcedures,     only : fatalError
  use particle_class,        only : particle, phaseCoord
  use particleDungeon_class, only : particleDungeon
  use outputFile_class,      only : outputFile

  use scoreMemory_class,     only : scoreMemory
  use tallyResult_class,     only : tallyResult, tallyResultEmpty
  use tallyClerk_inter,      only : tallyClerk

  implicit none
  private

  !!
  !! Simplest possible analog k-eff estimator that determines
  !! criticality by comparing population weight before and after transport cycle
  !! Assumes that constant to normalise fission site production is uniform in geometry
  !!
  !! SAMPLE DICTIOANRY INPUT:
  !!
  !! myClerk {
  !!   type keffAnalogClerk;
  !! }
  !!
  type, public,extends(tallyClerk) :: keffAnalogClerk
    private
    real(defReal) :: startPopWgt = ZERO ! Start of Generation Weight
    real(defReal) :: endPopWgt   = ZERO ! End of generation Weight
  contains
    ! Procedures used during build
    procedure :: init
    procedure :: validReports
    procedure :: getSize

    ! File reports and check status -> run-time procedures
    procedure :: reportCycleStart
    procedure :: reportCycleEnd

    ! Output procedures
    procedure  :: display
    procedure  :: print
    procedure  :: getResult
  end type keffAnalogClerk

  !!
  !! k-eff result class
  !!
  type,public, extends(tallyResult) :: keffResult
    real(defReal), dimension(2) :: keff = [ONE, ZERO]
  end type keffResult


contains

  !!
  !! Initialise from dictionary and name
  !!
  subroutine init(self, dict, name)
    class(keffAnalogClerk), intent(inout) :: self
    class(dictionary), intent(in)        :: dict
    character(nameLen), intent(in)       :: name
    character(100),parameter :: Here = 'init (keffAnalogClerk.f90)'

    ! Needs no settings, just load name
    call self % setName(name)

    ! Ensure correct initialisation to default values
    self % startPopWgt = ZERO
    self % endPopWgt   = ZERO

  end subroutine init

  !!
  !! Returns array of codes that represent diffrent reports
  !!
  function validReports(self) result(validCodes)
    class(keffAnalogClerk),intent(in)           :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    validCodes = [ cycleStart_CODE, cycleEnd_CODE ]

  end function validReports

  !!
  !! Return memory size of the clerk
  !!
  elemental function getSize(self) result(S)
    class(keffAnalogClerk), intent(in) :: self
    integer(shortInt)                  :: S

    S = 1

  end function getSize

  !!
  !! Process beginning of a cycle
  !!
  subroutine reportCycleStart(self, start, mem)
    class(keffAnalogClerk), intent(inout) :: self
    class(particleDungeon), intent(in)    :: start
    type(scoreMemory), intent(inout)      :: mem

    ! Update start population weight
    self % startPopWgt = self % startPopWgt + start % popWeight()


  end subroutine reportCycleStart

  !!
  !! Process end of the cycle
  !!
  subroutine reportCycleEnd(self, end, mem)
    class(keffAnalogClerk), intent(inout) :: self
    class(particleDungeon), intent(in)    :: end
    type(scoreMemory), intent(inout)      :: mem
    real(defReal)                         :: k_norm, k_eff

    ! Update end population weight
    self % endPopWgt = self % endPopWgt + end % popWeight()

    ! Close batch
    if( mem % lastCycle() ) then
      k_norm = end % k_eff

      ! Calculate and score analog estimate of k-eff
      k_eff =  self % endPopWgt / self % startPopWgt * k_norm
      call mem % accumulate(k_eff, self % getMemAddress() )

      self % startPopWgt = ZERO
      self % endPopWgt   = ZERO
    end if

  end subroutine reportCycleEnd

  !!
  !! Display convergance progress on the console
  !!
  subroutine display(self, mem)
    class(keffAnalogClerk), intent(in)  :: self
    type(scoreMemory), intent(in)       :: mem
    real(defReal)                       :: k, STD

    call mem % getResult(k, STD, self % getMemAddress())

    ! Print estimates to a console
    print '(A,F8.5,A,F8.5)', 'k-eff (analog): ',  k, ' +/- ', STD

  end subroutine display

  !!
  !! Write contents of the clerk in the slot to output file
  !!
  subroutine print(self, outFile, mem)
    class(keffAnalogClerk), intent(in) :: self
    class(outputFile), intent(inout)   :: outFile
    type(scoreMemory), intent(in)      :: mem
    real(defReal)                      :: k, STD
    character(nameLen)                 :: name

    ! Get result value
    call mem % getResult(k, STD, self % getMemAddress())

    ! Print to output file
    call outFile % startBlock(self % getName() )
    name = 'k_analog'
    call outFile % printResult(k, STD, name)
    call outFile % endBlock()

  end subroutine print

  !!
  !! Return result for interaction with Physics Package
  !! from the clerk in the slot
  !!
  pure subroutine getResult(self, res, mem)
    class(keffAnalogClerk), intent(in)              :: self
    class(tallyResult), allocatable, intent(inout)  :: res
    type(scoreMemory), intent(inout)                :: mem
    real(defReal)                                   :: k, STD

    ! Get result value
    call mem % getResult(k, STD, self % getMemAddress())

    allocate(res, source = keffResult([k, STD]))

  end subroutine getResult
    
end module keffAnalogClerk_class
