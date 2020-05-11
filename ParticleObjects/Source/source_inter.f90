module source_inter

  use numPrecision
  use particle_class,        only: particleState
  use particleDungeon_class, only: particleDungeon
  use dictionary_class,      only: dictionary
  use RNG_class,             only: RNG
  use geometry_inter,        only: geometry

  implicit none
  private

  !!
  !! Abstract interface of source for particles
  !!
  !! Source generates particles from specified distributions
  !! for, e.g., fixed source calcs. or to generate initial
  !! distribution for eigenvalue calcs
  !!
  !! In the simplest cases, all sampling procedures can be done
  !! by sampling from delta distributions, but this can be generalised
  !! to more complex distributions
  !!
  !! Private members:
  !!   geom -> Pointer to the geometry to ensure source is inside and
  !!           for more complicated source distribution
  !!
  !! Interface:
  !!   init              -> initialise the source
  !!   generate          -> generate particles to fill a dungeon
  !!   sampleParticle    -> sample particles from the corresponding distributions
  !!   sampleType        -> sets the particle type
  !!   samplePosition    -> samples the particle's position in the geometry
  !!   sampleEnergy      -> samples the particle's energy
  !!   sampleEnergyAngle -> samples the particle's energy and angle from corresponding distr.
  !!   kill              -> clean up the source
  !!
  type, public,abstract :: source
    private
    class(geometry), pointer, public       :: geom => null()
  contains
    procedure, non_overridable             :: generate
    procedure, non_overridable             :: sampleParticle
    procedure(init), deferred              :: init
    procedure(sampleType), deferred        :: sampleType
    procedure(samplePosition), deferred    :: samplePosition
    procedure(sampleEnergy), deferred      :: sampleEnergy
    procedure(sampleEnergyAngle), deferred :: sampleEnergyAngle
    procedure(kill), deferred              :: kill
  end type source

  abstract interface

    !!
    !! Initialise source from dictionary
    !!
    subroutine init(self, dict, geom)
      import :: source, &
                dictionary, &
                geometry, &
      class(source), intent(inout)         :: self
      class(dictionary), intent(in)        :: dict
      class(geometry), pointer, intent(in) :: geom
    end subroutine init

    !!
    !! sample particle type
    !!
    subroutine sampleType(self, p, rand)
      import :: source, &
                particleState, &
                RNG
      class(source), intent(inout)        :: self
      class(particleState), intent(inout) :: p
      class(RNG), pointer, intent(in)     :: rand
    end subroutine sampleType

    !!
    !! sample particle position
    !!
    subroutine samplePosition(self, p, rand)
      import :: source, &
                particleState, &
                RNG
      class(source), intent(inout)        :: self
      class(particleState), intent(inout) :: p
      class(RNG), pointer, intent(in)     :: rand
    end subroutine samplePosition

    !!
    !! sample particle energy
    !!
    subroutine sampleEnergy(self, p, rand)
      import :: source, &
                particleState, &
                RNG
      class(source), intent(inout)        :: self
      class(particleState), intent(inout) :: p
      class(RNG), pointer, intent(in)     :: rand
    end subroutine sampleEnergy
    
    !!
    !! sample particle energy-angle
    !!
    subroutine sampleEnergyAngle(self, p, rand)
      import :: source, &
                particleState, &
                RNG
      class(source), intent(inout)        :: self
      class(particleState), intent(inout) :: p
      class(RNG), pointer, intent(in)     :: rand
    end subroutine sampleEnergyAngle
    
    !!
    !! Deallocate memory used by source
    !!
    subroutine kill(self)
      import :: source
      class(source), intent(inout) :: self
    end subroutine kill

  end interface

contains

    !!
    !! Generate particles to populate a particleDungeon
    !!
    !! Fills a particle dungeon with n particles, sampled
    !! from the corresponding source distributions
    !!
    !! Args:
    !!   dungeon [inout] -> particle dungeon to be populated
    !!   n [in]          -> number of particles to place in dungeon
    !!
    !! Result:
    !!   A dungeon populated with n particles sampled from the source
    !!
    subroutine generate(self, dungeon, n, rand)
      class(source), intent(inout)         :: self
      type(particleDungeon), intent(inout) :: dungeon
      integer(shortInt), intent(in)        :: n
      class(RNG), pointer, intent(in)      :: rand
      type(particleState)                  :: p
      integer(shortInt)                    :: i

      ! Set dungeon size to begin
      call dungeon % setSize(n)

      ! Generate n particles to populate dungeon
      do i = 1, n
        p % wgt = ONE
        p % time = ZERO
        call self % sampleParticle(p, rand)
        call dungeon % replace(p, i)
      end do

    end subroutine generate

    !!
    !! Sample particle's phase space co-ordinates
    !!
    !! For the given source type, proceed through each distribution
    !! and sample the particle properties
    !!
    !! Args:
    !!   p [inout] -> particle to be over-written
    !!
    !! Result:
    !!   A particle sampled the prescribed source
    !!
    !! Errors:
    !!   Errors may occur in substituent sampling procedures
    !!
    subroutine sampleParticle(self, p, rand)
      class(source), intent(inout)       :: self
      type(particleState), intent(inout) :: p
      class(RNG), pointer, intent(in)    :: rand

      call self % sampleType(p, rand)
      call self % samplePosition(p, rand)
      call self % sampleEnergyAngle(p, rand)
      call self % sampleEnergy(p, rand)

    end subroutine sampleParticle

end module source_inter
