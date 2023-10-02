module source_inter

  use numPrecision
  use particle_class,        only : particleState
  use particleDungeon_class, only : particleDungeon
  use dictionary_class,      only : dictionary
  use RNG_class,             only : RNG
  use geometry_inter,        only : geometry

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
  !! Public members:
  !!   geom -> Pointer to the geometry to ensure source is inside and
  !!           for more complicated source distribution
  !!
  !! Interface:
  !!   init              -> initialise the source
  !!   generate          -> generate particles to fill a dungeon
  !!   sampleParticle    -> sample particles from the corresponding distributions
  !!   kill              -> clean up the source
  !!
  type, public,abstract :: source
    private
    class(geometry), pointer, public       :: geom => null()
  contains
    procedure, non_overridable             :: generate
    procedure(sampleParticle), deferred    :: sampleParticle
    procedure(init), deferred              :: init
    procedure(kill), deferred              :: kill
  end type source

  abstract interface

    !!
    !! Initialise source from dictionary & geometry
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
    !! Sample particle's phase space co-ordinates
    !!
    !! Generates a phase-space state for a single particle
    !!
    !! Args:
    !!   p [inout] -> particle to be over-written
    !!
    !! Result:
    !!   A particle sampled the prescribed source
    !!
    function sampleParticle(self, rand) result(p)
      import :: source, particleState, RNG
      class(source), intent(inout)       :: self
      class(RNG), intent(inout)          :: rand
      type(particleState)                :: p
    end function sampleParticle

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
      class(RNG), intent(in)               :: rand
      type(RNG), save                      :: pRand
      integer(shortInt)                    :: i
      !$omp threadprivate(pRand)

      ! Set dungeon size to begin
      call dungeon % setSize(n)

      ! Generate n particles to populate dungeon
      ! TODO: advance the rand after source generation!
      !       This should prevent reusing RNs during transport
      !$omp parallel
      pRand = rand
      !$omp do
      do i = 1, n
        call pRand % stride(i)
        call dungeon % replace(self % sampleParticle(pRand), i)
      end do
      !$omp end do
      !$omp end parallel

    end subroutine generate

    !!
    !! Return to uninitialised state
    !!
    elemental subroutine kill(self)
      class(source), intent(inout) :: self

      self % geom => null()

    end subroutine kill

end module source_inter
