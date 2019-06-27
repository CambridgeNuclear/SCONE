module centreOfMassClerk_class

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
  !!      type comClerk;
  !!      cycles 900;
  !!  }
  !!
  type, public, extends(tallyClerk) :: centreOfMassClerk
    private
    real(defReal),dimension(:,:),allocatable :: value            !! cycle-wise COM value
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
  end type centreOfMassClerk

contains

  !!
  !! Initialise clerk from dictionary and name
  !!
  subroutine init(self, dict, name)
    class(centreOfMassClerk), intent(inout)   :: self
    class(dictionary), intent(in)             :: dict
    character(nameLen), intent(in)            :: name

    ! Assign name
    call self % setName(name)

    ! Read number of cycles for which to track COM
    call dict % get(self % maxCycles, 'cycles')

    ! Allocate space for storing value
    allocate(self % value(self % maxCycles, 3))
    self % value = ZERO

  end subroutine init

  !!
  !! Returns array of codes that represent diffrent reports
  !!
  function validReports(self) result(validCodes)
    class(centreOfMassClerk),intent(in)        :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    validCodes = [cycleEnd_Code]

  end function validReports

  !!
  !! Return memory size of the clerk
  !!
  elemental function getSize(self) result(S)
    class(centreOfMassClerk), intent(in)   :: self
    integer(shortInt)                      :: S

    S = 3*self % maxCycles

  end function getSize

  !!
  !! Process cycle end
  !!
  subroutine reportCycleEnd(self, end, mem)
    class(centreOfMassClerk), intent(inout)   :: self
    class(particleDungeon), intent(in)        :: end
    type(scoreMemory), intent(inout)          :: mem
    type(particleState)                       :: state
    integer(shortInt)                         :: i, cc

    if ((self % currentCycle) < (self % maxCycles)) then

      self % currentCycle = self % currentCycle + 1
      cc = self % currentCycle

      ! Loop through population, scoring probabilities
      do i = 1,end % popSize()
        associate( state => end % get(i) )
          self % value(cc,:) = self % value(cc,:) + state % wgt * state % r
        end associate
      end do

      self % value(cc,:) = self % value(cc,:) / end % popWeight()

    end if

  end subroutine reportCycleEnd

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
    integer(shortInt)                    :: i
    real(defReal)                        :: val
    character(nameLen)                   :: name

    ! Begin block
    call outFile % startBlock(self % getName())

    ! Print COM
    name = 'CoMx'
    call outFile % startArray(name, [self % maxCycles])
    do i=1,self % maxCycles
      call outFile % addValue(self % value(i,1))
    end do
    call outFile % endArray()

    name = 'CoMy'
    call outFile % startArray(name, [self % maxCycles])
    do i=1,self % maxCycles
      call outFile % addValue(self % value(i,2))
    end do
    call outFile % endArray()

    name = 'CoMz'
    call outFile % startArray(name, [self % maxCycles])
    do i=1,self % maxCycles
      call outFile % addValue(self % value(i,3))
    end do
    call outFile % endArray()

    call outFile % endBlock()

  end subroutine print

  !!
  !! Returns to uninitialised state
  !!
  elemental subroutine kill(self)
    class(centreOfMassClerk), intent(inout) :: self

    if(allocated(self % value)) deallocate(self % value)
    self % currentCycle = 0
    self % maxCycles = 0

  end subroutine kill

end module centreOfMassClerk_class
