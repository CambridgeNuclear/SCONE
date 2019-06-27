module shannonEntropyClerk_class

  use numPrecision
  use tallyCodes
  use genericProcedures,          only : fatalError
  use dictionary_class,           only : dictionary
  use particle_class,             only : particle, particleState
  use particleDungeon_class,      only : particleDungeon
  use outputFile_class,           only : outputFile

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
  !! This is a prototype implementation of an implicit Shannon entropy estimator
  !! Takes cycle end reports to generate cycle-wise entropy
  !! Contains only a single map for discretisation
  !! Scores Shannon entropy for a user-specified number of cycles
  !!
  !! Notes:
  !!    ->
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
    class(tallyMap), allocatable             :: map
    real(defReal),dimension(:),allocatable   :: prob             !! probability of being in a given bin
    real(defReal),dimension(:),allocatable   :: value            !! cycle-wise value of entropy
    integer(shortInt)                        :: N = 0            !! Number of bins
    integer(shortInt)                        :: maxCycles = 0    !! Number of tally cycles
    integer(shortInt)                        :: currentCycle = 0 !! track current cycle

  contains
    ! Procedures used during build
    procedure  :: init
    procedure  :: validReports
    procedure  :: getSize

    ! File reports and check status -> run-time procedures
    procedure  :: reportCycleEnd

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

    ! Allocate space for storing entropy
    allocate(self % value(self % maxCycles))
    self % value = ZERO

    ! Read size of the map
    self % N = self % map % bins(0)

    ! Allocate space for storing probabilities
    allocate(self % prob(self % N))
    self % prob = ZERO

  end subroutine init

  !!
  !! Returns array of codes that represent diffrent reports
  !!
  function validReports(self) result(validCodes)
    class(shannonEntropyClerk),intent(in)      :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    validCodes = [cycleEnd_Code]

  end function validReports

  !!
  !! Return memory size of the clerk
  !!
  elemental function getSize(self) result(S)
    class(shannonEntropyClerk), intent(in) :: self
    integer(shortInt)                      :: S

    S = self % N + self % maxCycles

  end function getSize

  !!
  !! Process cycle end
  !!
  subroutine reportCycleEnd(self, end, mem)
    class(shannonEntropyClerk), intent(inout) :: self
    class(particleDungeon), intent(in)        :: end
    type(scoreMemory), intent(inout)          :: mem
    type(particleState)                       :: state
    integer(shortInt)                         :: i, j, cc, idx
    real(defReal)                             :: totWgt, one_log2

    if (self % currentCycle < self % maxCycles) then

      self % currentCycle = self % currentCycle + 1
      cc = self % currentCycle

      ! Loop through population, scoring probabilities
      do i = 1,end % popSize()
        associate( state => end % get(i) )
          idx = self % map % map(state)
          if( idx > 0) self % prob(idx) = self % prob(idx) + state % wgt
        end associate
      end do

      totWgt = end % popWeight()
      one_log2 = ONE/log(TWO)

      ! Loop through bins, summing entropy
      do j = 1,self % N
        self % prob(j) = self % prob(j)/totWgt
        if ((self % prob(j) > ZERO) .AND. (self % prob(j) < ONE)) then
          self % value(cc) = self % value(cc) - self % prob(j) * log(self % prob(j)) * one_log2
        end if
        self % prob(j) = ZERO
      end do

    end if

  end subroutine reportCycleEnd

  !!
  !! Display convergance progress on the console
  !!
  subroutine display(self, mem)
    class(shannonEntropyClerk), intent(in) :: self
    type(scoreMemory), intent(in)    :: mem

    print *, 'shannonEntropyClerk does not support display yet'

  end subroutine display

  !!
  !! Write contents of the clerk to output file
  !!
  subroutine print(self, outFile, mem)
    class(shannonEntropyClerk), intent(in) :: self
    class(outputFile), intent(inout) :: outFile
    type(scoreMemory), intent(in)    :: mem
    integer(shortInt)                :: i
    character(nameLen)               :: name

    ! Begin block
    call outFile % startBlock(self % getName())

    ! Print entropy
    name = 'shannonEntropy'

    call outFile % startArray(name, [self % maxCycles])

    do i=1,self % maxCycles
      call outFile % addValue(self % value(i))
    end do
    call outFile % endArray()

    call outFile % endBlock()

  end subroutine print

  !!
  !! Returns to uninitialised state
  !!
  elemental subroutine kill(self)
    class(shannonEntropyClerk), intent(inout) :: self

    if(allocated(self % map)) deallocate(self % map)
    if(allocated(self % prob)) deallocate(self % prob)
    if(allocated(self % value)) deallocate(self % value)
    self % N = 0
    self % currentCycle = 0
    self % maxCycles = 0

  end subroutine kill

end module shannonEntropyClerk_class
