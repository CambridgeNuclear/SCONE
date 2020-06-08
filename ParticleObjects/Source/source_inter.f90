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
  !! Extendable scource class procedures
  !!
  public :: kill 

  !!
  !! Abstract interface of source for particles
  !!
  !! Source generates particles from specified distributions
  !! for, e.g., fixed source calcs. or to generate initial
  !! distribution for eigenvalue calcs
  !!
  !! A single sample of a particle state is created by calling `sample****` procedures in order
  !! given in `sampleParticle` function. Note that `sampleEnergyAngle` is called AFTER
  !! `sampleEnergy`
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
    !! Args:
    !!   dict [in] -> dict containing point source information
    !!   geom [in] -> pointer to a geometry
    !!
    subroutine init(self, dict, geom)
      import :: source, &
                dictionary, &
                geometry
      class(source), intent(inout)         :: self
      class(dictionary), intent(in)        :: dict
      class(geometry), pointer, intent(in) :: geom
    end subroutine init

    !!
    !! Sample Type of a particle
    !!
    !! Sets 'Type' in the particle p (e.g. P_NEUTRON)
    !!
    !! Inputs:
    !!   p [inout] -> particle to be given a type
    !!   rand [in] -> random number generator
    !!
    subroutine sampleType(self, p, rand)
      import :: source, &
                particleState, &
                RNG
      class(source), intent(inout)        :: self
      class(particleState), intent(inout) :: p
      class(RNG), intent(inout)           :: rand
    end subroutine sampleType

    !!
    !! Sample particle position
    !!
    !! Sets position of the particle p
    !!
    !! Inputs:
    !!   p [inout] -> particle to be given a position
    !!   rand [in] -> random number generator
    !!
    subroutine samplePosition(self, p, rand)
      import :: source, &
                particleState, &
                RNG
      class(source), intent(inout)        :: self
      class(particleState), intent(inout) :: p
      class(RNG), intent(inout)           :: rand
    end subroutine samplePosition

    !!
    !! Sample particle Energy/Group
    !!
    !! Sets energy of a particle to a CE value or a MG index
    !! Also sets 'isMG' flag to .true. or .false.
    !!
    !! Inputs:
    !!   p [inout] -> particle to be given a position
    !!   rand [in] -> random number generator
    !!
    subroutine sampleEnergy(self, p, rand)
      import :: source, &
                particleState, &
                RNG
      class(source), intent(inout)        :: self
      class(particleState), intent(inout) :: p
      class(RNG), intent(inout)           :: rand
    end subroutine sampleEnergy

    !!
    !! Sample particle Energy/Group and angle Angle
    !!
    !! Sets diraction of a particle together with its energy.
    !! Sampling of energy is optional if Angle & Energy are uncorrelated
    !! Is called after `sampleEnergy`, to overwrite value provided by that subroutine
    !!
    !! Inputs:
    !!   p [inout] -> particle to be given a position
    !!   rand [in] -> random number generator
    !!
    subroutine sampleEnergyAngle(self, p, rand)
      import :: source, &
                particleState, &
                RNG
      class(source), intent(inout)        :: self
      class(particleState), intent(inout) :: p
      class(RNG), intent(inout)           :: rand
    end subroutine sampleEnergyAngle

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

    !!
    !! Return to uninitialised state
    !!
    subroutine kill(self)
      class(source), intent(inout) :: self

      self % geom => null()

    end subroutine kill

end module source_inter
