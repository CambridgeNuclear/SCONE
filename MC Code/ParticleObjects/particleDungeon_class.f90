module particleDungeon_class

  use numPrecision
  use genericProcedures,     only : fatalError, numToChar
  use particle_class,        only : particle, phaseCoord
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
    real(defReal)                 :: popWgt      ! Current population weight
    type(phaseCoord), dimension(:), allocatable :: prisoners  ! Data for the stored particles

  contains
    procedure  :: init
    procedure  :: isEmpty
    generic    :: detain  => detain_particle, detain_phaseCoord
    procedure  :: release
    procedure  :: normWeight
    procedure  :: normSize
    procedure  :: popSize
    procedure  :: popWeight
    procedure  :: printSourceToFile

    ! Private procedures
    procedure, private :: detain_particle
    procedure, private :: detain_phaseCoord
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
    self % popWgt = ZERO

  end subroutine init

  !!
  !! Store particle in the dungeon
  !!
  subroutine detain_particle(self,p)
    class(particleDungeon), intent(inout) :: self
    type(particle), intent(in)            :: p
    character(100),parameter              :: Here = 'detain_particle (particleDungeon_class.f90)'

    ! Increase population and weight
    self % pop = self % pop +1
    !self % popWgt = self % popWgt + p % w

    ! Check for population overflow
    if (self % pop > size(self % prisoners)) then
      print *, self % pop
      print *, size(self % prisoners)
      call fatalError(Here,'Run out of space for particles')
    end if

    ! Load new particle
    self % prisoners(self % pop) = p

  end subroutine detain_particle

  !!
  !! Store phaseCoord in the dungeon
  !!
  subroutine detain_phaseCoord(self,p_phase)
    class(particleDungeon), intent(inout) :: self
    type(phaseCoord), intent(in)          :: p_phase
    character(100), parameter             :: Here = 'detain_phaseCoord (particleDungeon_class.f90)'

    ! Increase population
    self % pop = self % pop +1
    !self % popWgt = self % popWgt + p_phase % wgt

    ! Check for population overflow
    if (self % pop > size(self % prisoners)) then
      print *, self % pop
      print *, size(self % prisoners)
      call fatalError(Here,'Run out of space for particles')
    end if

    ! Load new particle
    self % prisoners(self % pop) = p_phase

  end subroutine detain_phaseCoord

  !!
  !! Returns .true. if dungeon is empty
  !!
  function isEmpty(self) result(isIt)
    class(particleDungeon), intent(in) :: self
    logical(defBool)                   :: isIt

    isIt = (self % pop == 0)

  end function isEmpty


  !!
  !! Obtain the particle from the dungeon.
  !!
  subroutine release(self,p )
    class(particleDungeon), intent(inout) :: self
    type(particle), intent(inout)         :: p
    integer(shortInt)                     :: pop

    ! Load data into the particle
    pop = self % pop

    p = self % prisoners(pop)
    p % isDead = .false.

    ! Decrease population and weight
    self % pop = self % pop - 1
    !self % popWgt = self % popWgt - p % w
  end subroutine release

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
  !!
  subroutine normSize(self,N,rand)
    class(particleDungeon), intent(inout) :: self
    integer(shortInt), intent(in)         :: N
    class(RNG), intent(inout)             :: rand
    integer(shortInt)                     :: excessP
    integer(shortInt)                     :: i, idx

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

    ! Recalculate weight
    !self % popWgt = sum( self % prisoners(1:self % pop) % wgt )

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
  !! Returns total population weight
  !!
  function popWeight(self) result(wgt)
    class(particleDungeon), intent(in) :: self
    real(defReal)                      :: wgt

    wgt = sum( self % prisoners(1:self % pop) % wgt )
  end function popWeight

  !!
  !! Prints the position of fission sites to a file
  !! Used initially for looking at clustering
  !!
  subroutine printSourceToFile(self, name)
    class(particleDungeon), intent(in) :: self
    character(*), intent(in)           :: name
    character(256)                     :: filename
    integer(shortInt)                  :: i

    filename = trim(name)//'.txt'
    open(unit = 10, file = filename, status = 'new')

    ! Print out each particle co-ordinate
    do i = 1, self % pop
      write(10,'(A,A,A)') numToChar(self % prisoners(i) % r)
    end do

    ! Close the file
    close(10)

  end subroutine printSourceToFile

end module particleDungeon_class
