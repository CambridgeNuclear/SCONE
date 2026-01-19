module shannonEntropyClerk_class

  use numPrecision
  use tallyCodes
  use genericProcedures,          only : fatalError
  use display_func,               only : statusMsg
  use dictionary_class,           only : dictionary
  use particle_class,             only : particle, particleState
  use particleDungeon_class,      only : particleDungeon
  use outputFile_class,           only : outputFile
  use mpi_func,                   only : getMPIWorldSize

  ! Basic tally modules
  use scoreMemory_class,          only : scoreMemory
  use tallyClerk_inter,           only : tallyClerk

  ! Tally Maps
  use tallyMap_inter,             only : tallyMap
  use tallyMapFactory_func,       only : new_tallyMap

  implicit none
  private

  !!
  !! Shannon entropy estimator
  !!
  !! This is a prototype implementation of an implicit Shannon entropy estimator
  !! Takes cycle end reports to generate cycle-wise entropy
  !! Contains only a single map for discretisation
  !!
  !! Scores Shannon entropy for a user-specified number of cycles
  !!
  !! Sample dictionary input:
  !!
  !!  clerkName {
  !!      type shannonEntropyClerk;
  !!      map { <TallyMapDef> };
  !!      cycles 900;
  !!  }
  !!
  type, public, extends(tallyClerk) :: shannonEntropyClerk
    private
    !! Map defining the discretisation
    class(tallyMap), allocatable  :: map
    integer(shortInt)             :: N = 0            !! Number of bins
    integer(shortInt)             :: maxCycles = 0    !! Number of tally cycles
    integer(shortInt)             :: currentCycle = 0 !! track current cycle

  contains
    ! Procedures used during build
    procedure  :: init
    procedure  :: validReports
    procedure  :: getSize

    ! File reports and check status -> run-time procedures
    procedure  :: reportCycleEnd
    procedure  :: closeCycle

    ! Output procedures
    procedure  :: display
    procedure  :: print

    ! Deconstructor
    procedure  :: kill

  end type shannonEntropyClerk

contains

  !!
  !! Initialise clerk from dictionary and name
  !!
  subroutine init(self, dict, name)
    class(shannonEntropyClerk), intent(inout) :: self
    class(dictionary), intent(in)             :: dict
    character(nameLen), intent(in)            :: name

    ! Assign name
    call self % setName(name)

    ! Read map
    call new_tallyMap(self % map, dict % getDictPtr('map'))

    ! Read number of cycles for which to track entropy
    call dict % get(self % maxCycles, 'cycles')

    ! Read size of the map
    self % N = self % map % bins(0)

  end subroutine init

  !!
  !! Returns array of codes that represent diffrent reports
  !!
  function validReports(self) result(validCodes)
    class(shannonEntropyClerk),intent(in)      :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    validCodes = [cycleEnd_Code, closeCycle_CODE]

  end function validReports

  !!
  !! Return memory size of the clerk
  !!
  elemental function getSize(self) result(S)
    class(shannonEntropyClerk), intent(in) :: self
    integer(shortInt)                      :: S

    S = self % N + 1 + self % maxCycles

  end function getSize

  !!
  !! Process cycle end
  !!
  subroutine reportCycleEnd(self, end, mem)
    class(shannonEntropyClerk), intent(inout) :: self
    class(particleDungeon), intent(in)        :: end
    type(scoreMemory), intent(inout)          :: mem
    integer(shortInt)                         :: i, idx

    ! Increment cycle number
    self % currentCycle = self % currentCycle + 1

    if (self % currentCycle <= self % maxCycles) then

      ! Sample population weight for MPI
      call mem % score(end % popWeight(), self % getMemAddress())

      ! Loop through population, scoring probabilities
      do i = 1, end % popSize()

        associate(state => end % get(i))
          idx = self % map % map(state)
          if (idx == 0) cycle
          call mem % score(state % wgt, self % getMemAddress() + idx)
        end associate

      end do

    end if

  end subroutine reportCycleEnd

  !!
  !! Close cycle
  !!
  subroutine closeCycle(self, end, mem)
    class(shannonEntropyClerk), intent(inout) :: self
    class(particleDungeon), intent(in)        :: end
    type(scoreMemory), intent(inout)          :: mem
    integer(shortInt)                         :: i
    integer(longInt)                          :: ccIdx
    real(defReal)                             :: totWgt, one_log2, prob, val

    if (self % currentCycle <= self % maxCycles) then

      ! Get total population weight for this cycle
      totWgt = mem % getScore(self % getMemAddress())

      ! Initialise contants and counters
      val = ZERO
      one_log2 = ONE/log(TWO)

      ! Loop through bins, summing entropy
      do i = 1, self % N

        prob = mem % getScore(self % getMemAddress() + i)
        prob = prob / totWgt

        if ((prob > ZERO) .and. (prob < ONE)) then
          val = val - prob * log(prob) * one_log2
        end if

      end do

      ccIdx = self % getMemAddress() + self % N + self % currentCycle
      call mem % accumulate(val, ccIdx)

      ! Reset memory bins in preparation for next cycle
      do i = 0, self % N
        call mem % resetBin(self % getMemAddress() + i)
      end do

    end if

  end subroutine closeCycle

  !!
  !! Display convergance progress on the console
  !!
  subroutine display(self, mem)
    class(shannonEntropyClerk), intent(in) :: self
    type(scoreMemory), intent(in)    :: mem

    call statusMsg('shannonEntropyClerk does not support display yet')

  end subroutine display

  !!
  !! Write contents of the clerk to output file
  !!
  subroutine print(self, outFile, mem)
    class(shannonEntropyClerk), intent(in) :: self
    class(outputFile), intent(inout)       :: outFile
    type(scoreMemory), intent(in)          :: mem
    integer(shortInt)                      :: i, batches
    character(nameLen)                     :: name
    integer(longInt)                       :: ccIdx
    real(defReal)                          :: val

    ! Determine number of batches for normalisation with MPI
    if (mem % reduced) then
      batches = 1
    else
      batches = getMPIWorldSize()
    end if

    ! Begin block
    call outFile % startBlock(self % getName())

    ! Print entropy
    name = 'shannonEntropy'

    call outFile % startArray(name, [self % maxCycles])

    do i = 1, self % maxCycles
      ccIdx = self % getMemAddress() + self % N + i
      call mem % getResult(val, ccIdx, samples = batches)
      call outFile % addValue(val)
    end do

    call outFile % endArray()

    call outFile % endBlock()

  end subroutine print

  !!
  !! Returns to uninitialised state
  !!
  elemental subroutine kill(self)
    class(shannonEntropyClerk), intent(inout) :: self

    if (allocated(self % map)) deallocate(self % map)
    self % N = 0
    self % currentCycle = 0
    self % maxCycles = 0

  end subroutine kill

end module shannonEntropyClerk_class
