module eventClerk_class

  use numPrecision
  use tallyCodes
  use endfConstants
  use universalVariables
  use genericProcedures,          only : fatalError, numToChar
  use dictionary_class,           only : dictionary
  use particle_class,             only : particle, particleState
  use outputFile_class,           only : outputFile
  use scoreMemory_class,          only : scoreMemory
  use tallyClerk_inter,           only : tallyClerk, kill_super => kill

  ! Nuclear Data interface
  use nuclearDatabase_inter,      only : nuclearDatabase

  ! Tally Filters
  use tallyFilter_inter,          only : tallyFilter
  use tallyFilterFactory_func,    only : new_tallyFilter

  ! Tally Maps
  use tallyMap_inter,             only : tallyMap
  use tallyMapFactory_func,       only : new_tallyMap

  ! OpenMP functions
  use openmp_func,                only : ompGetMaxThreads, ompGetThreadNum

  implicit none
  private

  ! Memory size of the event object
  integer(shortInt), parameter :: eventSize = 12

  !!
  !! An object that records the state a particle was in during a
  !! particular collision.
  !!
  !! Stores data in reduced precision. Can be changed in future if necessary.
  !!
  type, private :: event
    ! Particle info
    real(defFlt), dimension(3) :: r
    real(defFlt), dimension(3) :: dir
    real(defFlt)               :: Eincident
    real(defFlt)               :: time
    real(defFlt)               :: w
    integer(shortInt)          :: brood
    
    ! Collision info
    real(defFlt)                :: Edeposit
    integer(shortInt)           :: MT
    !integer(shortInt)           :: ZAID ! Don't yet have available data in interface
  contains
    procedure :: build

  end type event

  !!
  !! Not a standard clerk: outputs particle events without averaging.
  !! No flux weighting is used. Events are stored locally, rather than
  !! in scoreMemory due to the unique 'score' type.
  !!
  !! Events are accumulated and output to a file when they reach a threshold
  !! number. These events contain information about a particle's state.
  !! They are written to a user-specified output file.
  !!
  !! Maps do not produce multiple files: they serve as filters deciding what 
  !! is output to the one chosen file.
  !!
  !! As this is used for tracking analog behaviour of neutrons, handleVirtual is
  !! off by default.
  !!
  !! The user can choose the write frequency and maximum number of events to write.
  !! This allows some control over computational expense and memory.
  !!
  !! TODO: allow filtering based on nuclide. May require interface change.
  !! TODO: allow maxEvents to be input directly (dictionary longInt restrictions)
  !!
  !! Private Members:
  !!   filter   -> Space to store tally Filter
  !!   map      -> Space to store tally Map
  !!   handleVirtual -> score on virtual collisions (probably not desirable)
  !!   outputFile -> path to the file where events are recorded
  !!   maxEvents -> maximum number of events; scoring stops when reached. 10M by default
  !!   frequency -> how many events before there is a file write/update. 500k by default
  !!   scores    -> array of events, with thread-privacy to avoid critical statements
  !!   lastEntry -> the last element of the scores array (for a given thread) which
  !!                was written to.
  !!   totalEvents -> total number of events which have been recorded
  !!
  !! Interface
  !!   tallyClerk Interface
  !!
  !! SAMPLE DICTIOANRY INPUT:
  !!
  !! myEventClerk {
  !!   type eventClerk;
  !!   file /home/sconeOutputs/myevents.txt;
  !!   maxScale 10;  ! Scales the maximum number of events since dictionary can't do longInt!
  !!   freq 5000000;
  !!   # handleVirtual 0; # default is 0   
  !!   # filter { <tallyFilter definition> } #
  !!   # map    { <tallyMap definition>    } #
  !! }
  !!
  type, public, extends(tallyClerk) :: eventClerk
    private
    ! Filter, Map & Vector of Responses
    class(tallyFilter), allocatable :: filter
    class(tallyMap), allocatable    :: map

    ! Settings
    character(pathLen) :: outputFile
    integer(longInt)   :: maxEvents = 100000000
    integer(shortInt)  :: frequency = 500000

    ! Storage
    type(event), dimension(:,:), allocatable     :: scores
    integer(shortInt), dimension(:), allocatable :: lastEntry
    integer(longInt)                             :: totalEvents = 0

    ! Settings
    logical(defBool)   :: handleVirtual = .false.

  contains
    ! Procedures used during build
    procedure  :: init
    procedure  :: kill
    procedure  :: validReports
    procedure  :: getSize

    ! File reports and check status -> run-time procedures
    procedure  :: reportOutColl

    ! Output procedures
    procedure           :: display
    procedure           :: print
    procedure, private  :: outputScores

  end type eventClerk

contains

  !!
  !! Initialise clerk from dictionary and name
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine init(self, dict, name)
    class(eventClerk), intent(inout) :: self
    class(dictionary), intent(in)    :: dict
    character(nameLen), intent(in)   :: name
    character(pathLen)               :: filepath
    integer(shortInt)                :: nThreads
    real(defReal)                    :: scaleMax
    character(100), parameter :: Here = 'init (eventClerk_class.f90)'

    ! Assign name
    call self % setName(name)

    ! Get output path
    call dict % get(filepath, 'file')
    self % outputFile = trim(filepath)

    call dict % getOrDefault(scaleMax, 'maxScale', ONE)
    self % maxEvents = int(self % maxEvents * scaleMax, longInt)
    if (self % maxEvents < 1) call fatalError(Here,&
            'The maximum number of events must be positive. Scale given: '&
            //numToChar(scaleMax))

    call dict % getOrDefault(self % frequency, 'freq', 5000000)
    if (self % frequency < 1) call fatalError(Here,&
            'The write frequency must be positive. Number given: '&
            //numToChar(self % frequency))

    ! Load filetr
    if( dict % isPresent('filter')) then
      call new_tallyFilter(self % filter, dict % getDictPtr('filter'))
    end if

    ! Load map
    if( dict % isPresent('map')) then
      call new_tallyMap(self % map, dict % getDictPtr('map'))
    end if

    ! Handle virtual collisions
    call dict % getOrDefault(self % handleVirtual,'handleVirtual', .false.)

    ! Size the scores array based on number of threads
    nThreads = ompGetMaxThreads()
    allocate(self % scores(self % frequency, nThreads))
    allocate(self % lastEntry(nThreads))
    self % lastEntry = 0

  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(eventClerk), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Kill and deallocate filter
    if (allocated(self % filter)) then
      deallocate(self % filter)
    end if

    ! Kill and deallocate map
    if (allocated(self % map)) then
      call self % map % kill()
      deallocate(self % map)
    end if

    self % handleVirtual = .false.
    self % outputFile = ''
    self % frequency = 5000000
    self % maxEvents = 100000000

  end subroutine kill

  !!
  !! Returns array of codes that represent different reports
  !!
  !! See tallyClerk_inter for details
  !!
  function validReports(self) result(validCodes)
    class(eventClerk),intent(in)               :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    validCodes = [outColl_CODE]

  end function validReports

  !!
  !! Return memory size of the clerk
  !!
  !! See tallyClerk_inter for details
  !!
  elemental function getSize(self) result(S)
    class(eventClerk), intent(in) :: self
    integer(shortInt)                 :: S

    S = 0

  end function getSize

  !!
  !! Process outgoing collision report
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportOutColl(self, p, MT, muL, xsData, mem)
    class(eventClerk), intent(inout)      :: self
    class(particle), intent(in)           :: p
    integer(shortInt), intent(in)         :: MT
    real(defReal), intent(in)             :: muL
    class(nuclearDatabase), intent(inout) :: xsData
    type(scoreMemory), intent(inout)      :: mem
    type(particleState)                   :: state
    integer(shortInt)                     :: binIdx, bin, id
    character(100), parameter :: Here = 'reportOutColl (eventClerk_class.f90)'

    if (self % totalEvents >= self % maxEvents) return
    if (MT == noInteraction .and. .not. self % handleVirtual) return

    ! Get current particle state
    state = p

    ! Check if within filter
    if (allocated(self % filter)) then
      if (self % filter % isFail(state)) return
    end if

    ! Find bin index
    if (allocated(self % map)) then
      binIdx = self % map % map(state)
    else
      binIdx = 1
    end if

    ! Return if invalid bin index
    if (binIdx == 0) return

    ! Convert the particle to an event
    ! Store the event in the thread-private bank
    id = ompGetThreadNum() + 1
    bin = self % lastEntry(id) + 1
    call self % scores(bin, id) % build(p, MT)
    self % lastEntry(id) = bin

    ! If the bank is full, flush it to output
    if (bin == self % frequency) then
      !$omp critical
      ! Escape in case all events were completed while waiting on
      ! the critical
      if (self % totalEvents >= self % maxEvents) then
        self % lastEntry(id) = 0
      else

        call self % outputScores(id)
      
        ! Increment the total number of events
        self % totalEvents = self % totalEvents + self % lastEntry(id)
        self % lastEntry(id) = 0
      end if
      !$omp end critical
    end if

  end subroutine reportOutColl

  !!
  !! Display convergence progress on the console
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine display(self, mem)
    class(eventClerk), intent(in) :: self
    type(scoreMemory), intent(in) :: mem

    print *, 'eventClerk does not support display yet'

  end subroutine display

  !!
  !! Write contents of the clerk to output file.
  !! In this case, only the map is written to the usual output.
  !! Meanwhile, the events are written to the event file if they
  !! haven't already been flushed.
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine print(self, outFile, mem)
    class(eventClerk), intent(in)              :: self
    class(outputFile), intent(inout)           :: outFile
    type(scoreMemory), intent(in)              :: mem
    integer(shortInt)                          :: i

    ! Make sure that score arrays have been flushed already.
    ! If not, write the remainder to output.
    do i = 1, size(self % lastEntry)
      if (self % lastEntry(i) > 0) call self % outputScores(i)
    end do

    ! Write out any map that was used to the usual output file,
    ! but scores are written out separately.

    ! Begin block
    call outFile % startBlock(self % getName())

    ! If collision clerk has map print map information
    if (allocated(self % map)) then
      call self % map % print(outFile)
    end if
    
    call outFile % endBlock()

  end subroutine print

  !!
  !! Send events in a thread's bank to the output file
  !!
  subroutine outputScores(self, id)
    class(eventClerk), intent(in) :: self
    integer(shortInt), intent(in) :: id
    integer(shortInt)             :: i, u, ios

    associate(ev => self % scores(:,id))


      ! Open the file and append the event info
      open(newunit=u, file=trim(self % outputFile), status='unknown', &
           position='append', action='write', iostat=ios)

      do i = 1, self % lastEntry(id)
        ! Consider changing to a binary output in future
        ! Would require a corresponding python script to translate
        write(u,'(*(g0,1x))') ev(i) % r, ev(i) % dir, ev(i) % Eincident, &
                ev(i) % time, ev(i) % w, ev(i) % brood, ev(i) % Edeposit, &
                ev(i) % MT
      end do
      close(u)

    end associate

  end subroutine outputScores

!!
!! Event procedures
!!

  !!
  !! Set an event given a particle and an MT
  !!
  subroutine build(self, p, MT)
    class(event), intent(inout)   :: self
    type(particle), intent(in)    :: p
    integer(shortInt), intent(in) :: MT

    self % r = real(p % rGlobal(), defFlt)
    self % dir = real(p % dirGlobal(), defFlt)
    self % Eincident = real(p % preCollision % E, defFlt)
    self % time = real(p % time, defFlt)
    self % w = real(p % w, defFlt)
    self % brood = p % broodID
    self % MT = MT

    ! Particle was absorbed
    if (p % isDead) then
      self % Edeposit = real(p % precollision % E, defFlt)

    ! Scatter occurred - compute energy change from pre-collision info
    else
      self % Edeposit = real(p % preCollision % E - p % E, defFlt)
    end if
    
  end subroutine build

end module eventClerk_class
