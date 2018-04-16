module particleDungeon_class

  use numPrecision
  use genericProcedures,     only : fatalError
  use particle_class,        only : particle
  use RNG_class,             only : RNG

  implicit none
  private

  !!
  !! Support type to contain all data relevant to the particle detained in the dungeon
  !!
  type, private :: phaseCoord
    real(defReal)              :: wgt  = 0.0     ! Particle weight
    real(defReal),dimension(3) :: r    = 0.0     ! Global position
    real(defReal),dimension(3) :: dir  = 0.0     ! Global direction
    real(defReal)              :: E    = 0.0     ! Energy
    integer(shortInt)          :: G    = 0       ! Energy group
    logical(defBool)           :: isMG = .false. ! Is neutron multi-group
  contains
    generic  :: assignment(=) => phaseCoord_fromParticle
    procedure, private :: phaseCoord_fromParticle
  end type


  !!
  !! particleDungeon stores particle phase-space
  !! Used in eigenvalue calculation to store fission sites generated in a cycle
  !! Similar structures are refered to as:
  !! Store: MONK and Serpent(?)
  !! Fission Bank: OpenMC and MCNP(?)
  !!
  !! For convinience it allows to store value of k-eff that can be retrieved to adjust fission site
  !! generation rate during a calculation.
  !!
  !! Technically not the whole particles are stored but only the key data defined in phaseCoord
  !!
  !! Dungeon works like stack:
  !! throw(particle)   -> adda a particle to the top
  !! release(particle) -> removes a particle from the top
  !! isEmpty()         -> returns .true. if there are no more particles
  !! init(maxSize)     -> allocate space to store maximum of maxSize particles
  !! normWeight(totWgt)-> normalise dungeon population so its total weight is totWgt
  !! normSize(N)       -> normalise dungeon population so it contains N particles
  !!
  type, public :: particleDungeon
    private
    real(defReal),public          :: k_eff = 1.0 ! k-eff for fission site generation rate normalisation
    integer(shortInt)             :: pop         ! Current population size of the dungeon
    type(phaseCoord), dimension(:), allocatable :: prisoners  ! Data for the stored particles

  contains
    procedure  :: init
    procedure  :: isEmpty
    procedure  :: throw
    procedure  :: release
    procedure  :: normWeight
    procedure  :: normSize
    procedure  :: popSize

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
  !! Throw particle into the dungeon (store it)
  !!
  subroutine throw(self,p)
    class(particleDungeon), intent(inout) :: self
    type(particle), intent(in)            :: p
    character(100),parameter              :: Here = 'throw (particleDungeon_class.f90)'

    ! Increase population
    self % pop = self % pop +1

    ! Check for population overflow
    if (self % pop > size(self % prisoners)) then
      call fatalError(Here,'Run out of space for particles')
    end if

    ! Load new particle
    self % prisoners(self % pop) = p

  end subroutine throw

  !!
  !! Returns .true. if dungeon is empty
  !!
  function isEmpty(self) result(isIt)
    class(particleDungeon), intent(in) :: self
    logical(defBool)                   :: isIt

    isIt = (self % pop == 0)

  end function isEmpty


  !!
  !! Obtain the particle from the dungeon. **** Will need to be change for more rebust loading
  !!
  subroutine release(self,p )
    class(particleDungeon), intent(inout) :: self
    type(particle), intent(inout)         :: p
    integer(shortInt)                     :: pop

    ! Load data into the particle
    pop = self % pop
    call p % teleport(self % prisoners(pop) % r)
    call p % point(self % prisoners(pop) % dir)
    p % w    = self % prisoners(pop) % wgt
    p % E    = self % prisoners(pop) % E
    p % G    = self % prisoners(pop) % G
    p % isMg = self % prisoners(pop) % isMG
    p % isDead = .false.

    ! Decrease population
    self % pop = self % pop - 1

  end subroutine release

  !!
  !! Normalise total weight of the particles in the dungeon to match provided value
  !!
  subroutine normWeight(self,totWgt)
    class(particleDungeon), intent(inout) :: self
    real(defReal), intent(in)             :: totWgt
    real(defReal)                         :: factor

    ! Behold the glory of Fortran! Normalisation of weights in two lines
    factor = totWgt / sum(self % prisoners % wgt)
    self % prisoners % wgt = self % prisoners % wgt * factor

  end subroutine normWeight

  !!
  !! Normalise total number of particles in the dungeon to match the provided number
  !! Randomly duplicate or remove particles to match the number
  !!
  subroutine normSize(self,N,rand)
    class(particleDungeon), intent(inout) :: self
    integer(shortInt), intent(in)         :: N
    class(RNG), intent(inout)             :: rand
    integer(shortInt)                     :: excessP
    integer(shortInt)                     :: i, idx
    real(defReal)                         :: r1
    real(defReal)                         :: flipProb

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
  !! Returns number of neutrons in the dungeon
  !!
  function popSize(self) result(pop)
    class(particleDungeon), intent(in) :: self
    integer(shortInt)                  :: pop

    pop = self % pop

  end function popSize

  !!
  !! Copy particle into phase coordinates
  !!
  subroutine phaseCoord_fromParticle(LHS,RHS)
    class(phaseCoord), intent(out)  :: LHS
    type(particle), intent(in)      :: RHS

    LHS % wgt  = RHS % w
    LHS % r    = RHS % rGlobal()
    LHS % dir  = RHS % dirGlobal()
    LHS % E    = RHS % E
    LHS % G    = RHS % G
    LHS % isMG = RHS % isMG

  end subroutine phaseCoord_fromParticle




end module particleDungeon_class
