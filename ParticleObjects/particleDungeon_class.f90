module particleDungeon_class

  use numPrecision
  use genericProcedures,     only : fatalError, numToChar, swap
  use particle_class,        only : particle, particleState
  use RNG_class,             only : RNG
  use heapQueue_class,       only : heapQueue

  use mpi_func,              only : isMPIMaster, getMPIWorldSize, getMPIRank, getOffset
#ifdef MPI
  use mpi_func,              only : mpi_gather, mpi_allgather, mpi_send, mpi_recv, &
                                    mpi_Bcast, MPI_COMM_WORLD, MASTER_RANK, MPI_DEFREAL, &
                                    MPI_SHORTINT, MPI_LONGINT, MPI_PARTICLE_STATE, &
                                    MPI_STATUS_IGNORE, particleStateDummy
#endif

  implicit none
  private

  !!
  !! particleDungeon stores particle phase-space
  !! Used in eigenvalue calculation to store fission sites generated in a cycle
  !! Similar structures are refered to as:
  !! Store: MONK and Serpent(?)
  !! Fission Bank: OpenMC and MCNP(?)
  !!
  !! NOTE INCONSISTENT DEFINITIONS
  !! ****
  !! For convenience it allows storing the value of k-eff that can be retrieved to adjust fission site
  !! generation rate during a calculation. It is not currently used during normalisation but
  !! is used by analog k-eff tally. It is necessary to clarify behaviour.
  !! ****
  !!
  !! Technically not the whole particles are stored but only the key data defined in phaseCoord
  !!
  !! Dungeon can work like stacks or arrays. Stack-like behaviour is not really thread safe
  !! so it can be utilised when collecting and processing secondary particles in history
  !! that should be processed during the course of one cycle. Alternatively, one can use the
  !! critical variations of the stack-like procedures.
  !! Array-like behaviour allows to easily distribute particles among threads. As long as indices
  !! assigned to different threads do not overlap, reading is thread-safe (I hope-MAK).
  !!
  !!
  !! INTERFACE:
  !!   Stack-like interface:
  !!     detain(particle)          -> add a particle to the top
  !!     detainCritical(particle)  -> add a particle to the top with a critical operation
  !!     release(particle)         -> removes a particle from the top. Sets p % isDead = .false.
  !!     releaseCritical(particle) -> removes a particle from the top with a critical operation.
  !!                                   Sets p % isDead = .false.
  !!
  !!   Array-like interface:
  !!     replace(particle, i) -> overwrite prisoner data at index i
  !!     copy(particle, i)    -> copy prisoner at index i into particle. Sets p % isDead = .false.
  !!     get(i)               -> function returns particle state at index i
  !!
  !!   Misc procedures:
  !!     isEmpty()         -> returns .true. if there are no more particles
  !!     cleanPop()        -> kill or prisoners
  !!     normWeight(totWgt)-> normalise dungeon population so its total weight is totWgt
  !!     normSize(N)       -> normalise dungeon population so it contains N particles
  !!                          does not take nonuniform weight of particles into account
  !!     setSize(n)        -> sizes dungeon to have n dummy particles for ease of overwriting
  !!     printToFile(name) -> prints population in ASCII format to file "name"
  !!
  !!   Build procedures:
  !!     init(maxSize)     -> allocate space to store maximum of maxSize particles
  !!     kill()            -> return to uninitialised state
  !!
  type, public :: particleDungeon
    private
    real(defReal), public :: k_eff = ONE ! k-eff for fission site generation rate normalisation
    integer(shortInt)     :: pop = 0     ! Current population size of the dungeon

    ! Storage space
    type(particleState), dimension(:), allocatable, public :: prisoners

  contains
    !! Build procedures
    procedure  :: init
    procedure  :: kill

    !! Stack-like interface
    generic    :: detain => detain_particle, detain_particleState
    generic    :: detainCritical => detainCritical_particle, detainCritical_particleState
    procedure  :: release
    procedure  :: releaseCritical

    !! Array-like interface
    generic    :: replace => replace_particle, replace_particleState
    procedure  :: copy
    procedure  :: get

    !! Misc Procedures
    procedure  :: isEmpty
    procedure  :: normWeight
    procedure  :: normSize
    procedure  :: cleanPop
    procedure  :: popSize
    procedure  :: popWeight
    procedure  :: setSize
    procedure  :: printToFile
    procedure  :: sortByBroodID
    procedure  :: samplingWoReplacement

    ! Private procedures
    procedure, private :: detain_particle
    procedure, private :: detain_particleState
    procedure, private :: detainCritical_particle
    procedure, private :: detainCritical_particleState
    procedure, private :: replace_particle
    procedure, private :: replace_particleState
#ifdef MPI
    procedure, private :: loadBalancing
#endif

  end type particleDungeon

contains

  !!
  !! Allocate space for the particles
  !!
  subroutine init(self,maxSize)
    class(particleDungeon), intent(inout) :: self
    integer(shortInt), intent(in)         :: maxSize

    if (allocated(self % prisoners)) deallocate(self % prisoners)
    allocate(self % prisoners(maxSize))
    self % pop    = 0

  end subroutine init

  !!
  !! Deallocate memory and return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(particleDungeon), intent(inout) :: self

    ! Reset settings
    self % pop = 0

    ! Deallocate memeory
    if (allocated(self % prisoners)) deallocate(self % prisoners)

  end subroutine kill

  !!
  !! Store particle in the dungeon
  !!
  subroutine detain_particle(self,p)
    class(particleDungeon), intent(inout) :: self
    class(particle), intent(in)           :: p
    integer(shortInt)                     :: pop
    character(100),parameter              :: Here = 'detain_particle (particleDungeon_class.f90)'

    !$omp atomic capture
    ! Increase population and weight
    self % pop = self % pop + 1
    pop = self % pop
    !$omp end atomic

    ! Check for population overflow
    if (pop > size(self % prisoners)) then
      call fatalError(Here,'Run out of space for particles.&
                           & Max size:'//numToChar(size(self % prisoners)) //&
                            ' Current population: ' // numToChar(self % pop))
    end if

    ! Load new particle
    self % prisoners(pop) = p

  end subroutine detain_particle

  !!
  !! Store particle in the dungeon with a critical operation
  !!
  subroutine detainCritical_particle(self,p)
    class(particleDungeon), intent(inout) :: self
    class(particle), intent(in)           :: p
    integer(shortInt)                     :: pop
    character(100),parameter              :: Here = 'detainCritical_particle (particleDungeon_class.f90)'

    !$omp critical (dungeon)
    ! Increase population and weight
    self % pop = self % pop + 1
    pop = self % pop

    ! Check for population overflow
    if (pop > size(self % prisoners)) then
      call fatalError(Here,'Run out of space for particles.&
                           & Max size:'//numToChar(size(self % prisoners)) //&
                            ' Current population: ' // numToChar(self % pop))
    end if

    ! Load new particle
    self % prisoners(pop) = p
    !$omp end critical (dungeon)

  end subroutine detainCritical_particle

  !!
  !! Store phaseCoord in the dungeon
  !!
  subroutine detain_particleState(self,p_state)
    class(particleDungeon), intent(inout) :: self
    type(particleState), intent(in)       :: p_state
    integer(shortInt)                     :: pop
    character(100), parameter    :: Here = 'detain_particleState (particleDungeon_class.f90)'

    ! Increase population
    !$omp atomic capture
    self % pop = self % pop + 1
    pop = self % pop
    !$omp end atomic

    ! Check for population overflow
    if (pop > size(self % prisoners)) then
      call fatalError(Here,'Run out of space for particles.&
                           & Max size:'//numToChar(size(self % prisoners)) //&
                            ' Current population: ' // numToChar(self % pop))
    end if

    ! Load new particle
    self % prisoners(pop) = p_state

  end subroutine detain_particleState

  !!
  !! Store phaseCoord in the dungeon with a critical operation
  !!
  subroutine detainCritical_particleState(self,p_state)
    class(particleDungeon), intent(inout) :: self
    type(particleState), intent(in)       :: p_state
    integer(shortInt)                     :: pop
    character(100), parameter    :: Here = 'detainCritical_particleState (particleDungeon_class.f90)'

    ! Increase population
    !$omp critical (dungeon)
    self % pop = self % pop + 1
    pop = self % pop

    ! Check for population overflow
    if (pop > size(self % prisoners)) then
      call fatalError(Here,'Run out of space for particles.&
                           & Max size:'//numToChar(size(self % prisoners)) //&
                            ' Current population: ' // numToChar(self % pop))
    end if

    ! Load new particle
    self % prisoners(pop) = p_state
    !$omp end critical (dungeon)

  end subroutine detainCritical_particleState

  !!
  !! Pop the particle from the top of the dungeon.
  !! Makes particle alive at exit
  !!
  subroutine release(self, p)
    class(particleDungeon), intent(inout) :: self
    type(particle), intent(inout)         :: p
    integer(shortInt)                     :: pop

    !$omp atomic capture
    ! Decrease population
    pop = self % pop
    self % pop = self % pop - 1
    !$omp end atomic

    ! Load data into the particle
    p = self % prisoners(pop)
    p % isDead = .false.

  end subroutine release

  !!
  !! Pop the particle from the top of the dungeon with a critical operation.
  !! Makes particle alive at exit
  !!
  subroutine releaseCritical(self, p)
    class(particleDungeon), intent(inout) :: self
    type(particle), intent(inout)         :: p
    integer(shortInt)                     :: pop

    !$omp critical (dungeon)
    ! Decrease population
    pop = self % pop
    self % pop = self % pop - 1

    ! Load data into the particle
    p = self % prisoners(pop)
    !$omp end critical (dungeon)

    p % isDead = .false.

  end subroutine releaseCritical

  !!
  !! Replace data of particle prisoner at the index idx with particle
  !!
  subroutine replace_particle(self, p, idx)
    class(particleDungeon), intent(inout) :: self
    class(particle), intent(in)           :: p
    integer(shortInt), intent(in)         :: idx
    character(100),parameter :: Here = 'replace_particle (particleDungeon_class.f90)'

    ! Protect against out-of-bounds access
    if (idx <= 0 .or. idx > self % pop) then
      call fatalError(Here,'Out of bounds access with idx: '// numToChar(idx)// &
                           ' with particle population of: '// numToChar(self % pop))
    end if

    ! Load new particle
    self % prisoners(idx) = p

  end subroutine replace_particle

  !!
  !! Replace data of particle prisoner at the index idx with phaseCoords
  !!
  subroutine replace_particleState(self, p, idx)
    class(particleDungeon), intent(inout) :: self
    type(particleState), intent(in)       :: p
    integer(shortInt), intent(in)         :: idx
    character(100),parameter :: Here = 'replace_particleState (particleDungeon_class.f90)'

    ! Protect against out-of-bounds access
    if (idx <= 0 .or. idx > self % pop) then
      call fatalError(Here,'Out of bounds access with idx: '// numToChar(idx)// &
                           ' with particle population of: '// numToChar(self % pop))
    end if

    ! Load new particle
    self % prisoners(idx) = p

  end subroutine replace_particleState

  !!
  !! Copy particle from a location inside the dungeon
  !!
  !! Makes particle alive at exit. Also sets the broodID of the particle
  !! making it ready to be transported.
  !!
  !! Args:
  !!  p   [inout] -> Particle to be filled with data
  !!  idx [in]    -> Index of the particle to be copied
  !!
  !! Errors:
  !!  fatalError if requested index is 0, -ve or above current population
  !!
  subroutine copy(self, p, idx)
    class(particleDungeon), intent(in) :: self
    type(particle), intent(inout)      :: p
    integer(shortInt), intent(in)      :: idx
    character(100), parameter :: Here = 'copy (particleDungeon_class.f90)'

    ! Protect against out-of-bounds access
    if (idx <= 0 .or. idx > self % pop) then
      call fatalError(Here,'Out of bounds access with idx: '// numToChar(idx)// &
                           ' with particle population of: '// numToChar(self % pop))
    end if

    ! Load data into the particle
    p = self % prisoners(idx)
    p % isDead = .false.
    p % broodID = idx

  end subroutine copy

  !!
  !! Return particleState from a location inside the dungeon
  !! Gives fatalError if requested index is 0, -ve or above current population
  !!
  function get(self, idx) result(state)
    class(particleDungeon), intent(in) :: self
    integer(shortInt), intent(in)      :: idx
    type(particleState)                :: state
    character(100), parameter :: Here = 'get (particleDungeon_class.f90)'

    ! Protect against out-of-bounds access
    if( idx <= 0 .or. idx > self % pop ) then
      call fatalError(Here,'Out of bounds access with idx: '// numToChar(idx)// &
                           ' with particle population of: '// numToChar(self % pop))
    end if

    ! Explicit copy. Will be changed soon
    state = self % prisoners(idx)

  end function get

  !!
  !! Returns .true. if dungeon is empty
  !!
  function isEmpty(self) result(isIt)
    class(particleDungeon), intent(in) :: self
    logical(defBool)                   :: isIt

    isIt = (self % pop == 0)

  end function isEmpty

  !!
  !! Normalise total weight of the particles in the dungeon to match provided value
  !!
  subroutine normWeight(self,totWgt)
    class(particleDungeon), intent(inout) :: self
    real(defReal), intent(in)             :: totWgt
    real(defReal)                         :: factor

    ! Behold the glory of Fortran! Normalisation of weights in two lines
    factor = totWgt / sum(self % prisoners(1:self % pop) % wgt)
    self % prisoners % wgt = self % prisoners % wgt * factor

  end subroutine normWeight

  !!
  !! Normalise total number of particles in the dungeon to match the provided number
  !!
  !! Does not take weight of a particle into account!
  !!
  subroutine normSize(self, totPop, rand)
    class(particleDungeon), intent(inout) :: self
    integer(shortInt), intent(in)         :: totPop
    class(RNG), intent(inout)             :: rand
    type(RNG)                             :: rankRand, masterRand
    type(heapQueue)                       :: heap
    real(defReal)                         :: threshold, rn
    integer(longInt)                      :: seedTemp
    integer(shortInt)                     :: maxbroodID, totSites, excess, heapSize, &
                                             n_duplicates, i, j, n_copies, count, nRanks, &
                                             rank
    integer(longInt), dimension(:), allocatable  :: seeds
    integer(shortInt), dimension(:), allocatable :: keepers, popSizes
#ifdef MPI
    integer(shortInt)                            :: error
#endif
    character(100), parameter :: Here = 'normSize (particleDungeon_class.f90)'

    ! Determine the maximum brood ID and sort the dungeon for OMP reproducibility
    maxBroodID = maxval(self % prisoners(1:self % pop) % broodID)
    call self % sortByBroodID(maxbroodID)

    ! Get MPI world size and allocate rng seed vector, needed by all processes
    nRanks = getMPIWorldSize()
    allocate(seeds(nRanks), popSizes(nRanks))

    ! Initialise popSizes with the correct value for when only one process is used
    popSizes  = self % pop
    seeds     = 0
    threshold = ONE

#ifdef MPI
    ! Get the population sizes of all ranks into the array popSizes in master branch
    call mpi_gather(self % pop, 1, MPI_SHORTINT, popSizes, 1, MPI_SHORTINT, MASTER_RANK, MPI_COMM_WORLD, error)
#endif

    ! In the master process, calculate sampling threshold for the whole population
    ! and send it to all processes
    if (isMPIMaster()) then

      ! Calculate number of sites generated over all processes and difference to target population
      totSites = sum(popSizes)
      excess   = totSites - totPop

      ! Assign heapQueue size according to case, accounting for cases where the whole
      ! population might have to be replicated due to massive undersampling
      if (excess < 0) then
        n_duplicates = modulo(-excess, totSites)
        heapSize     = n_duplicates
      else
        heapSize = excess
      end if

      if (heapSize /= 0) then

        ! Copy rng
        masterRand = rand

        ! Initialise heapQueue and push upper bound larger than 1.0
        call heap % init(heapSize)
        call heap % pushReplace(TWO)

        ! Loop to generate totSites random numbers to fill the heapQueue
        do i = 1, nRanks

          ! Save rng seed: this will be the starting seed in the i-th rank
          seeds(i) = masterRand % currentState()

          ! Populate heapQueue
          do j = 1, popSizes(i)
            rn = masterRand % get()
            if (rn < heap % maxValue()) call heap % pushReplace(rn)
          end do

        end do

        ! Save sampling threshold
        threshold = heap % maxValue()

      end if

    end if

    ! Broadcast threshold, excess and random number seeds to all processes
#ifdef MPI
    call MPI_Bcast(threshold, 1, MPI_DEFREAL, MASTER_RANK, MPI_COMM_WORLD)
    call MPI_Bcast(excess, 1, MPI_SHORTINT, MASTER_RANK, MPI_COMM_WORLD)
    call MPI_Bcast(seeds, nRanks, MPI_LONGINT, MASTER_RANK, MPI_COMM_WORLD)
#endif

    ! Get local process rank and initialise local rng with the correct seed
    rank = getMPIRank()
    seedTemp = seeds(rank + 1)
    call rankRand % init(seedTemp)

    ! Perform the actual sampling
    if (excess > 0) then

      allocate(keepers(self % pop))
      keepers = 0
      count   = 0

      ! Loop over source sites
      do i = 1, self % pop

        ! Save the indexes of the sites to keep and increment counter
        if (rankRand % get() > threshold) then
          count = count + 1
          keepers(count) = i
        end if

      end do

      ! Loop through accepted sites to save them
      do i = 1, count
        self % prisoners(i) = self % prisoners(keepers(i))
      end do

      ! Update population number
      self % pop = count

    elseif (excess < 0) then

      ! Check if copies have to be made and the number of particles to duplicate with the sampling
      totSites     = excess + totPop
      n_copies     = -excess / totSites
      n_duplicates = modulo(-excess, totSites)

      ! Copy all the particles the number of times needed
      do i = 1, n_copies
        self % prisoners(self % pop * i + 1 : self % pop * (i + 1)) = self % prisoners(1:self % pop)
      end do

      ! Loop over population to duplicate from
      count = self % pop * (n_copies + 1)

      if (n_duplicates /= 0) then

        do i = 1, self % pop
          ! Save duplicated particles at the end of the dungeon
          if (rankRand % get() <= threshold) then
            count = count + 1
            self % prisoners(count) = self % prisoners(i)
          end if
        end do

      end if

      ! Update population number
      self % pop = count

      ! Determine the maximum brood ID and sort the dungeon again for MPI reproducibility
      maxBroodID = maxval(self % prisoners(1:self % pop) % broodID)
      call self % sortByBroodID(maxbroodID)

    end if

    ! Get the new population in the case of one thread
    popSizes = self % pop

#ifdef MPI
    ! Get the updated population numbers from all processes
    call mpi_allgather(self % pop, 1, MPI_SHORTINT, popSizes, 1, MPI_SHORTINT, MPI_COMM_WORLD, error)
#endif

    ! Check that normalisation worked
    if (sum(popSizes) /= totPop) call fatalError(Here, 'Normalisation failed!')

#ifdef MPI
    ! Perform load balancing by redistributing particles across processes
    if (nRanks > 1) call self % loadBalancing(totPop, nRanks, rank, popSizes)
#endif

  end subroutine normSize

  !!
  !! Perform nearest neighbor load balancing
  !!
#ifdef MPI
  subroutine loadBalancing(self, totPop, nRanks, rank, popSizes)
    class(particleDungeon), intent(inout)        :: self
    integer(shortInt), intent(in)                :: totPop
    integer(shortInt), intent(in)                :: nRanks
    integer(shortInt), intent(in)                :: rank
    integer(shortInt), dimension(:), intent(in)  :: popSizes
    integer(shortInt), dimension(:), allocatable :: rankOffsets
    integer(shortInt), dimension(2)              :: offset, targetOffset
    integer(shortInt)                            :: mpiOffset, excess, error
    type(particleState), dimension(:), allocatable      :: stateBuffer
    type(particleStateDummy), dimension(:), allocatable :: dummyBuffer

    ! Get expected particle population in each process via the offset
    mpiOffset = getOffset(totPop)

    ! Communicates the offset from all processes to all processes
    allocate(rankOffsets(nRanks))
    call mpi_allgather(mpiOffset, 1, MPI_SHORTINT, rankOffsets, 1, MPI_SHORTINT, MPI_COMM_WORLD, error)

    ! Calculate actual and target cumulative number of sites in the processes before
    offset(1)       = sum(popSizes(1 : rank))
    offset(2)       = sum(popSizes(1 : rank + 1))
    targetOffset(1) = rankOffsets(rank + 1)
    if (rank + 1 == nRanks) then
      targetOffset(2) = totPop
    else
      targetOffset(2) = rankOffsets(rank + 2)
    end if

    ! If needed, send/receive particle states from/to the end of the dungeon
    excess = offset(2) - targetOffset(2)

    if (excess > 0) then

      ! Send particles from the end of the dungeon to the rank above
      stateBuffer = self % prisoners(self % pop - excess + 1 : self % pop)
      call initStateDummies(stateBuffer, dummyBuffer)
      call mpi_send(dummyBuffer, excess, MPI_PARTICLE_STATE, rank + 1, rank, MPI_COMM_WORLD, error)
      self % pop = self % pop - excess

    elseif (excess < 0) then

      ! Receive particles from the rank above and store them at the end of the dungeon
      excess = -excess
      allocate(dummyBuffer(excess))
      call mpi_recv(dummyBuffer, excess, MPI_PARTICLE_STATE, rank + 1, rank + 1, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, error)
      call createStatesFromDummies(dummyBuffer, stateBuffer)
      self % prisoners(self % pop + 1 : self % pop + excess) = stateBuffer
      self % pop = self % pop + excess

    end if

    if (excess /= 0) deallocate(stateBuffer, dummyBuffer)

    ! If needed, send/receive particle states from/to the beginning of the dungeon
    excess = offset(1) - targetOffset(1)

    if (excess < 0) then

      ! Send particles from the beginning of the dungeon to the rank below
      excess = -excess
      stateBuffer = self % prisoners(1 : excess)
      call initStateDummies(stateBuffer, dummyBuffer)
      call mpi_send(dummyBuffer, excess, MPI_PARTICLE_STATE, rank - 1, rank, MPI_COMM_WORLD, error)

      ! Move the remaining particles to the beginning of the dungeon
      self % prisoners(1 : self % pop - excess) = self % prisoners(excess + 1 : self % pop)
      self % pop = self % pop - excess

    elseif (excess > 0) then

      ! Receive particles from the rank below and store them at the beginning of the dungeon
      allocate(dummyBuffer(excess))
      call mpi_recv(dummyBuffer, excess, MPI_PARTICLE_STATE, rank - 1, rank - 1, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, error)
      call createStatesFromDummies(dummyBuffer, stateBuffer)
      self % prisoners(excess + 1 : self % pop + excess) = self % prisoners(1 : self % pop)
      self % prisoners(1 : excess) = stateBuffer
      self % pop = self % pop + excess

    end if

  end subroutine loadBalancing
#endif

  !!
  !! Helper procedure for MPI loadBalancing
  !!
  !! Given an input array of particleState, it allocates a vector of particleStateDummy
  !! of the same length and copies all the particleState variables. It returns
  !! the array of particleStateDummy
  !!
#ifdef MPI
  subroutine initStateDummies(states, dummies)
    type(particleState), dimension(:), intent(in)                    :: states
    type(particleStateDummy), dimension(:), allocatable, intent(out) :: dummies
    integer(shortInt)                                                :: i, N

    ! Allocate particleStateDummy array
    N = size(states)
    allocate(dummies(N))

    ! Copy particleState attributes
    do i = 1, N
      dummies(i) % r   = states(i) % r
      dummies(i) % dir = states(i) % dir
    end do

    dummies % wgt  = states % wgt
    dummies % E    = states % E
    dummies % G    = states % G
    dummies % isMG = states % isMG
    dummies % type = states % type
    dummies % time = states % time
    dummies % matIdx     = states % matIdx
    dummies % uniqueID   = states % uniqueID
    dummies % cellIdx    = states % cellIdx
    dummies % collisionN = states % collisionN
    dummies % broodID    = states % broodID

  end subroutine initStateDummies
#endif

  !!
  !! Helper procedure for MPI loadBalancing
  !!
  !! Given an input array of particleStateDummy, it allocates a vector of particleState
  !! of the same length and copies all the particleStateDummy variables. It returns
  !! the array of particleState
  !!
#ifdef MPI
  subroutine createStatesFromDummies(dummies, states)
    type(particleStateDummy), dimension(:), intent(in)          :: dummies
    type(particleState), dimension(:), allocatable, intent(out) :: states
    integer(shortInt)                                           :: i, N

    ! Allocate particleState array
    N = size(dummies)
    allocate(states(N))

    ! Copy particleStateDummy attributes
    do i = 1, N
      states(i) % r   = dummies(i) % r
      states(i) % dir = dummies(i) % dir
    end do

    states % wgt  = dummies % wgt
    states % E    = dummies % E
    states % G    = dummies % G
    states % isMG = dummies % isMG
    states % type = dummies % type
    states % time = dummies % time
    states % matIdx     = dummies % matIdx
    states % uniqueID   = dummies % uniqueID
    states % cellIdx    = dummies % cellIdx
    states % collisionN = dummies % collisionN
    states % broodID    = dummies % broodID

  end subroutine createStatesFromDummies
#endif


  !!
  !! Normalise total number of particles in the dungeon to match the provided number.
  !! Randomly duplicate or remove particles to match the number, performing
  !! sampling without replacement (Knuth S algorithm).
  !!
  !! Does not take weight of a particle into account!
  !!
  !! NOTE: this procedure is not currently used
  !!
  subroutine samplingWoReplacement(self, N, rand)
    class(particleDungeon), intent(inout) :: self
    integer(shortInt), intent(in)         :: N
    class(RNG), intent(inout)             :: rand
    integer(shortInt)                     :: excessP, n_copies, n_duplicates
    integer(shortInt)                     :: i, idx, maxBroodID
    integer(shortInt), dimension(:), allocatable :: duplicates
    character(100), parameter :: Here = 'samplingWoReplacement (particleDungeon_class.f90)'

    ! Protect against invalid N
    if (N > size(self % prisoners)) then
      call fatalError(Here,'Requested size: '//numToChar(N) //&
                           'is greater then max size: '//numToChar(size(self % prisoners)))
    else if (N <= 0) then
      call fatalError(Here,'Requested size: '//numToChar(N) //' is not +ve')
    end if

    ! Determine the maximum brood ID and sort the dungeon
    maxBroodID = maxval(self % prisoners(1:self % pop) % broodID)
    call self % sortByBroodID(maxbroodID)

    ! Calculate excess particles to be removed
    excessP = self % pop - N

    if (excessP > 0) then ! Reduce population with reservoir sampling
      do i = N + 1, self % pop
        ! Select new index. Copy data if it is in the safe zone (<= N).
        idx = int(i * rand % get()) + 1
        if (idx <= N) then
          self % prisoners(idx) = self % prisoners(i)
        end if
      end do
      self % pop = N

    else if (excessP < 0) then ! Clone randomly selected particles
      ! For a massive undersampling duplicate (or n-plicate) particles
      excessP = -excessP
      n_copies = excessP / self % pop
      n_duplicates = modulo(excessP, self % pop)

      ! Copy all particles the maximum possible number of times
      do i = 1, n_copies
        self % prisoners(self % pop * i + 1 : self % pop * (i + 1)) = self % prisoners(1:self % pop)
      end do

      ! Choose the remainder particles to duplicate without replacement
      duplicates = [(i, i = 1, n_duplicates)]
      do i = n_duplicates + 1, self % pop
        idx = int(i * rand % get()) + 1
        if (idx <= n_duplicates) then
          duplicates(idx) = i
        end if
      end do
      self % pop = self % pop * (n_copies + 1)

      ! Copy the duplicated particles at the end
      do i = 1, n_duplicates
        self % prisoners(self % pop + i) = self % prisoners(duplicates(i))
      end do
      self % pop = N

    end if

  end subroutine samplingWoReplacement

  !!
  !! Reorder the dungeon so the brood ID is in the ascending order
  !!
  !! Args:
  !!   k [in] -> Maximum brood ID
  !!
  subroutine sortByBroodID(self, k)
    class(particleDungeon), intent(inout)        :: self
    integer(shortInt), intent(in)                :: k
    integer(shortInt), dimension(k)              :: count
    integer(shortInt)                            :: i, id, loc, c, j
    integer(shortInt), dimension(:), allocatable :: perm
    type(particleState)                          :: tmp
    character(100), parameter :: Here = 'sortBybroodID (particleDungeon_class.f90)'

    ! Count number of particles with each brood ID
    count = 0
    do i = 1, self % pop
      id = self % prisoners(i) % broodID
      if (id < 1 .or. id > k) call fatalError(Here, 'Brood ID out of range: '//numToChar(id))
      count(id) = count(id) + 1
    end do

    ! Convert to starting index
    loc = 1
    do i = 1, k
      c = count(i)
      count(i) = loc
      loc = loc + c
    end do

    ! Create the permutation array
    allocate(perm(self % pop))
    do i = 1, self % pop
      id = self % prisoners(i) % broodID
      loc = count(id)
      count(id) = count(id) + 1
      perm(loc) = i
    end do

    ! Permute particles
    do i = 1, self % pop
      j = i

      do while(i /= perm(i))
        loc = perm(i)

        ! Swap elements
        if (loc /= j) then
          tmp = self % prisoners(j)
          self % prisoners(j) = self % prisoners(loc)
          self % prisoners(loc) = tmp
        end if
        call swap(perm(i), perm(loc))

        j = loc

      end do

    end do

  end subroutine sortByBroodID


  !!
  !! Kill or particles in the dungeon
  !!
  pure subroutine cleanPop(self)
    class(particleDungeon), intent(inout) :: self

    self % pop = 0

  end subroutine cleanPop

  !!
  !! Returns number of neutrons in the dungeon
  !!
  function popSize(self) result(pop)
    class(particleDungeon), intent(in) :: self
    integer(shortInt)                  :: pop

    pop = self % pop

  end function popSize

  !!
  !! Returns total population weight
  !!
  function popWeight(self) result(wgt)
    class(particleDungeon), intent(in) :: self
    real(defReal)                      :: wgt

    wgt = sum(self % prisoners(1:self % pop) % wgt)

  end function popWeight

  !!
  !! Set size of the dungeon to n
  !!
  !! Sets population to arbitrary size n
  !! All stored particles revert to default initialisation state
  !!
  !! Args:
  !!   n [in] -> Requested size of the population
  !!
  !! Errors:
  !!   fatalError if n is invalid (not +ve)
  !!
  subroutine setSize(self, n)
    class(particleDungeon), intent(inout) :: self
    integer(shortInt), intent(in)         :: n
    character(100), parameter :: Here = 'setSize (particleDungeon_class.f90)'

    if (n <= 0) call fatalError(Here, 'Requested population is not +ve: '//numToChar(n))

    ! Set population
    self % pop = n

    ! Make sure enough space is avaliable
    if (allocated(self % prisoners)) then
      if (size(self % prisoners) < n) then
        deallocate(self % prisoners)
        allocate(self % prisoners(n))
      end if

    else
      allocate(self % prisoners(n))
    end if

    ! Set known (default) state to all particles
    call self % prisoners % kill()

  end subroutine setSize

  !!
  !! Prints the position of fission sites to a file
  !! Used initially for looking at clustering
  !!
  subroutine printToFile(self, name)
    class(particleDungeon), intent(in) :: self
    character(*), intent(in)           :: name
    character(256)                     :: filename
    integer(shortInt)                  :: i

    filename = trim(name)//'.txt'
    open(unit = 10, file = filename, status = 'new')

    ! Print out each particle co-ordinate
    do i = 1, self % pop
      write(10, *) self % prisoners(i) % r, self % prisoners(i) % dir, &
                   self % prisoners(i) % E, self % prisoners(i) % G, &
                   self % prisoners(i) % broodID
    end do

    ! Close the file
    close(10)

  end subroutine printToFile

end module particleDungeon_class
