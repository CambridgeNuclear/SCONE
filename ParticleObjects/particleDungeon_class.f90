module particleDungeon_class

  use numPrecision
  use genericProcedures,     only : fatalError, numToChar
  use particle_class,        only : particle, particleState
  use RNG_class,             only : RNG

  implicit none
  private

  !!
  !! particleDungeon stores particle phase-space
  !! Used in eigenvalue calculation to store fission sites generated in a cycle
  !! Similar structures are refered to as:
  !! Store: MONK and Serpent(?)
  !! Fission Bank: OpenMC and MCNP(?)
  !!
  !! NOTE INCONSISTANT DEFINITIONS
  !! ****
  !! For convinience it allows to store value of k-eff that can be retrieved to adjust fission site
  !! generation rate during a calculation. It is not currently used during normalisation but
  !! is used by analog k-eff tally. It is necessary to clarify behaviour.
  !! ****
  !!
  !! Technically not the whole particles are stored but only the key data defined in phaseCoord
  !!
  !! Dungeon can work like stacks or arrays. Stack-like behaviour is not really thread safe
  !! so it can be utilised when collecting and processing secondary particles in history
  !! that should be processed during the course of one cycle. Array-like behaviour allows to
  !! easily distribute particles among threads. As long as indices assign to diffrent threads
  !! do not overlap, reading is thread-safe (I hope-MAK).
  !!
  !!
  !! INTERFACE:
  !!   Stack-like interface:
  !!     detain(particle)   -> adda a particle to the top
  !!     release(particle)  -> removes a particle from the top. Sets p % isDead = .false.
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
  !!                          does not take ununiform weight of particles into account
  !!     setSize(n)        -> sizes dungeon to have n dummy particles for ease of overwriting
  !!     printToFile(name) -> prints population in ASCII format to file "name"
  !!
  !!   Build procedures:
  !!     init(maxSize)     -> allocate space to store maximum of maxSize particles
  !!     kill()            -> return to uninitialised state
  !!
  type, public :: particleDungeon
    private
    real(defReal),public :: k_eff = ONE ! k-eff for fission site generation rate normalisation
    integer(shortInt)    :: pop = 0     ! Current population size of the dungeon

    ! Storage space
    type(particleState), dimension(:), allocatable :: prisoners

  contains
    !! Build procedures
    procedure  :: init
    procedure  :: kill

    !! Stack-like interface
    generic    :: detain  => detain_particle, detain_particleState
    procedure  :: release

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

    ! Private procedures
    procedure, private :: detain_particle
    procedure, private :: detain_particleState
    procedure, private :: replace_particle
    procedure, private :: replace_particleState
  end type particleDungeon

contains

  !!
  !! Allocate space for the particles
  !!
  subroutine init(self,maxSize)
    class(particleDungeon), intent(inout) :: self
    integer(shortInt), intent(in)         :: maxSize

    if(allocated(self % prisoners)) deallocate(self % prisoners)
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
    if(allocated(self % prisoners)) deallocate(self % prisoners)

  end subroutine kill

  !!
  !! Store particle in the dungeon
  !!
  subroutine detain_particle(self,p)
    class(particleDungeon), intent(inout) :: self
    class(particle), intent(in)           :: p
    character(100),parameter              :: Here = 'detain_particle (particleDungeon_class.f90)'

    ! Increase population and weight
    self % pop = self % pop +1

    ! Check for population overflow
    if (self % pop > size(self % prisoners)) then
      call fatalError(Here,'Run out of space for particles.&
                           & Max size:'//numToChar(size(self % prisoners)) //&
                            ' Current population: ' // numToChar(self % pop))
    end if

    ! Load new particle
    self % prisoners(self % pop) = p

  end subroutine detain_particle

  !!
  !! Store phaseCoord in the dungeon
  !!
  subroutine detain_particleState(self,p_state)
    class(particleDungeon), intent(inout) :: self
    type(particleState), intent(in)       :: p_state
    character(100), parameter    :: Here = 'detain_particleState (particleDungeon_class.f90)'

    ! Increase population
    self % pop = self % pop +1

    ! Check for population overflow
    if (self % pop > size(self % prisoners)) then
      call fatalError(Here,'Run out of space for particles.&
                           & Max size:'//numToChar(size(self % prisoners)) //&
                            ' Current population: ' // numToChar(self % pop))
    end if

    ! Load new particle
    self % prisoners(self % pop) = p_state

  end subroutine detain_particleState

  !!
  !! Pop the particle from the top of the dungeon.
  !! Makes particle alive at exit
  !!
  subroutine release(self, p)
    class(particleDungeon), intent(inout) :: self
    type(particle), intent(inout)         :: p
    integer(shortInt)                     :: pop

    ! Load data into the particle
    pop = self % pop
    p = self % prisoners(pop)
    p % isDead = .false.

    ! Decrease population
    self % pop = self % pop - 1

  end subroutine release

  !!
  !! Replace data of particle prisoner at the index idx with particle
  !!
  subroutine replace_particle(self, p, idx)
    class(particleDungeon), intent(inout) :: self
    class(particle), intent(in)           :: p
    integer(shortInt), intent(in)         :: idx
    character(100),parameter :: Here = 'relplace_particle (particleDungeon_class.f90)'

    ! Protect agoinst out-of-bounds acces
    if( idx <= 0 .or. idx > self % pop ) then
      call fatalError(Here,'Out of bounds acces with idx: '// numToChar(idx)// &
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
    character(100),parameter :: Here = 'relplace_particleState (particleDungeon_class.f90)'

    ! Protect agoinst out-of-bounds acces
    if( idx <= 0 .or. idx > self % pop ) then
      call fatalError(Here,'Out of bounds acces with idx: '// numToChar(idx)// &
                           ' with particle population of: '// numToChar(self % pop))
    end if

    ! Load new particle
    self % prisoners(idx) = p

  end subroutine replace_particleState


  !!
  !! Copy particle from a location inside the dungeon
  !! Makes particle alive at exit
  !! Gives fatalError if requested index is 0, -ve or above current population
  !!
  subroutine copy(self, p, idx)
    class(particleDungeon), intent(in) :: self
    type(particle), intent(inout)      :: p
    integer(shortInt), intent(in)      :: idx
    character(100), parameter :: Here = 'copy (particleDungeon_class.f90)'

    ! Protect agoinst out-of-bounds acces
    if( idx <= 0 .or. idx > self % pop ) then
      call fatalError(Here,'Out of bounds acces with idx: '// numToChar(idx)// &
                           ' with particle population of: '// numToChar(self % pop))
    end if

    ! Load data into the particle
    p = self % prisoners(idx)
    p % isDead = .false.

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

    ! Protect agoinst out-of-bounds acces
    if( idx <= 0 .or. idx > self % pop ) then
      call fatalError(Here,'Out of bounds acces with idx: '// numToChar(idx)// &
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
  !! Randomly duplicate or remove particles to match the number
  !! Does not take weight of a particle into account!
  !!
  subroutine normSize(self,N,rand)
    class(particleDungeon), intent(inout) :: self
    integer(shortInt), intent(in)         :: N
    class(RNG), intent(inout)             :: rand
    integer(shortInt)                     :: excessP
    integer(shortInt)                     :: i, idx
    character(100), parameter :: Here =' normSize (particleDungeon_class.f90)'

    ! Protect against invalid N
    if( N > size(self % prisoners)) then
      call fatalError(Here,'Requested size: '//numToChar(N) //&
                           'is greather then max size: '//numToChar(size(self % prisoners)))
    else if ( N <= 0 ) then
      call fatalError(Here,'Requested size: '//numToChar(N) //' is not +ve')
    end if

    ! Calculate excess particles to be removed
    excessP = self % pop - N

    if (excessP > 0 ) then ! Reduce population with reservoir sampling
      do i=N,self % pop
        ! Select new index. Copy data if it is in the safe zone (<= N).
        idx = int(i * rand % get())+1
        if (idx <= N) then
          self % prisoners(idx) = self % prisoners(i)
        end if
      end do
      self % pop = N

    else if (excessP < 0) then ! Clone randomly selected particles
      do i = self % pop, N
        idx = int(self % pop * rand % get()) + 1
        self % prisoners(i) = self % prisoners(idx)
      end do
      self % pop = N

    end if

  end subroutine normSize

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

    wgt = sum( self % prisoners(1:self % pop) % wgt )

  end function popWeight

  !!
  !! Set population to have a particular size
  !! Used for ease of overwriting, e.g., when sampling from a source 
  !!
  subroutine setSize(self, n)
    class(particleDungeon), intent(inout) :: self
    integer(shortInt), intent(in)         :: n

    self % pop = n
    if (size(self % prisoners) < 3*n) then
      deallocate(self % prisoners)
      allocate(self % prisoners(3*n))
    end if

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
      write(10,'(8A)') numToChar(self % prisoners(i) % r), &
                       numToChar(self % prisoners(i) % dir), &
                       numToChar(self % prisoners(i) % E), &
                       numToChar(self % prisoners(i) % G)
    end do

    ! Close the file
    close(10)

  end subroutine printToFile

end module particleDungeon_class
