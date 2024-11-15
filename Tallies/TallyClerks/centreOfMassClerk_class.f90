module centreOfMassClerk_class

  use numPrecision
  use tallyCodes
  use genericProcedures,          only : fatalError, numToChar
  use dictionary_class,           only : dictionary
  use particle_class,             only : particle, particleState
  use particleDungeon_class,      only : particleDungeon
  use outputFile_class,           only : outputFile
  use mpi_func,                   only : getMPIWorldSize

  ! Basic tally modules
  use scoreMemory_class,          only : scoreMemory
  use tallyClerk_inter,           only : tallyClerk

  implicit none
  private

  !!
  !! Centre of Mass estimator
  !! This is a prototype implementation of an implicit neutron COM
  !! Takes cycle end reports to generate cycle-wise COM in three co-ordinates
  !! Contains only a single map for discretisation
  !!
  !! Notes:
  !!    ->dim1 of value is cycle
  !!    ->dim2 of value is dimension (x, y, z)
  !! Sample dictionary input:
  !!
  !!  clerkName {
  !!      type centreOfMassClerk;
  !!      cycles 900;
  !!  }
  !!
  type, public, extends(tallyClerk) :: centreOfMassClerk
    private
    integer(shortInt) :: maxCycles = 0    !! Number of tally cycles
    integer(shortInt) :: currentCycle = 0 !! track current cycle

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

  end type centreOfMassClerk

contains

  !!
  !! Initialise clerk from dictionary and name
  !!
  subroutine init(self, dict, name)
    class(centreOfMassClerk), intent(inout)   :: self
    class(dictionary), intent(in)             :: dict
    character(nameLen), intent(in)            :: name
    character(100), parameter :: Here = 'init (centreOfMassClerk_class.f90)'

    ! Assign name
    call self % setName(name)

    ! Read number of cycles for which to track COM
    call dict % get(self % maxCycles, 'cycles')

    ! Check on cycle number
    if (self % maxCycles <= 0) then
      call fatalError(Here, 'Number of cycles shuold be positive. It is: '//trim(numToChar(self % maxCycles)))
    end if

  end subroutine init

  !!
  !! Returns array of codes that represent diffrent reports
  !!
  function validReports(self) result(validCodes)
    class(centreOfMassClerk),intent(in)        :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    validCodes = [cycleEnd_Code, closeCycle_CODE]

  end function validReports

  !!
  !! Return memory size of the clerk
  !!
  elemental function getSize(self) result(S)
    class(centreOfMassClerk), intent(in)   :: self
    integer(shortInt)                      :: S

    S = 3 * self % maxCycles + 1

  end function getSize

  !!
  !! Process cycle end
  !!
  subroutine reportCycleEnd(self, end, mem)
    class(centreOfMassClerk), intent(inout)   :: self
    class(particleDungeon), intent(in)        :: end
    type(scoreMemory), intent(inout)          :: mem
    integer(shortInt)                         :: i
    integer(longInt)                          :: cc
    real(defReal), dimension(3)               :: val

    ! Increment cycle number
    self % currentCycle = self % currentCycle + 1

    if ((self % currentCycle) <= (self % maxCycles)) then

      cc = self % currentCycle

      ! Loop through population, scoring probabilities
      do i = 1, end % popSize()

        associate(state => end % get(i))
          val = state % wgt * state % r
          call mem % score(val(1), self % getMemAddress() + 3*(cc - 1))
          call mem % score(val(2), self % getMemAddress() + 3*(cc - 1) + 1)
          call mem % score(val(3), self % getMemAddress() + 3*(cc - 1) + 2)
        end associate

      end do

      ! Sample population weight for MPI
      call mem % score(end % popWeight(), self % getMemAddress() + 3*self % maxCycles)

    end if

  end subroutine reportCycleEnd

  !!
  !! Process cycle end
  !!
  subroutine closeCycle(self, end, mem)
    class(centreOfMassClerk), intent(inout)   :: self
    class(particleDungeon), intent(in)        :: end
    type(scoreMemory), intent(inout)          :: mem
    real(defReal)                             :: norm
    integer(longInt)                          :: cc

    if ((self % currentCycle) <= (self % maxCycles)) then

      ! Retrieve population weight and use it to normalise scores
      norm = mem % getScore(self % getMemAddress() + 3*self % maxCycles)
      norm = ONE / norm

      cc = self % currentCycle

      ! Make sure results don't get normalised arbitrarily
      call mem % closeBin(norm, self % getMemAddress() + 3*(cc - 1))
      call mem % closeBin(norm, self % getMemAddress() + 3*(cc - 1) + 1)
      call mem % closeBin(norm, self % getMemAddress() + 3*(cc - 1) + 2)

    end if

  end subroutine closeCycle

  !!
  !! Display convergance progress on the console
  !!
  subroutine display(self, mem)
    class(centreOfMassClerk), intent(in) :: self
    type(scoreMemory), intent(in)        :: mem

    print *, 'centreOfMassClerk does not support display yet'

  end subroutine display

  !!
  !! Write contents of the clerk to output file
  !!
  subroutine print(self, outFile, mem)
    class(centreOfMassClerk), intent(in) :: self
    class(outputFile), intent(inout)     :: outFile
    type(scoreMemory), intent(in)        :: mem
    integer(shortInt)                    :: i, batches
    character(nameLen)                   :: name
    integer(longInt)                     :: ccIdx
    real(defReal)                        :: val

    ! Determine number of batches for normalisation with MPI
    if (mem % reduced) then
      batches = 1
    else
      batches = getMPIWorldSize()
    end if

    ! Begin block
    call outFile % startBlock(self % getName())

    ! Print COM
    name = 'CoMx'
    call outFile % startArray(name, [self % maxCycles])
    do i = 1, self % maxCycles
      ccIdx = self % getMemAddress() + 3*(i - 1)
      call mem % getResult(val, ccIdx, samples = batches)
      call outFile % addValue(val)
    end do
    call outFile % endArray()

    name = 'CoMy'
    call outFile % startArray(name, [self % maxCycles])
    do i = 1, self % maxCycles
      ccIdx = self % getMemAddress() + 3*(i - 1) + 1
      call mem % getResult(val, ccIdx, samples = batches)
      call outFile % addValue(val)
    end do
    call outFile % endArray()

    name = 'CoMz'
    call outFile % startArray(name, [self % maxCycles])
    do i = 1, self % maxCycles
      ccIdx = self % getMemAddress() + 3*(i - 1) + 2
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
    class(centreOfMassClerk), intent(inout) :: self

    self % currentCycle = 0
    self % maxCycles = 0

  end subroutine kill

end module centreOfMassClerk_class
