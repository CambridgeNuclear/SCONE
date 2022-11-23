module source_inter

  use numPrecision
  use particle_class,        only : particle, particleState
  use particleDungeon_class, only : particleDungeon
  use dictionary_class,      only : dictionary
  use RNG_class,             only : RNG
  use geometry_inter,        only : geometry
  use genericProcedures,     only : fatalError

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
  !!   append            -> generate new particles to add to an existing dungeon
  !!   sampleParticle    -> sample particles from the corresponding distributions
  !!   kill              -> clean up the source
  !!
  type, public,abstract :: source
    private
    class(geometry), pointer, public       :: geom => null()
  contains
    procedure                              :: generate
    procedure                              :: append
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
      class(RNG), intent(inout)            :: rand
      integer(shortInt)                    :: i

      ! Set dungeon size to begin
      call dungeon % setSize(n)

      ! Generate n particles to populate dungeon
      do i = 1, n
        call dungeon % replace(self % sampleParticle(rand), i)
      end do

    end subroutine generate

    !!
    !! Generate particles to add to a particleDungeon without overriding
    !! particles already present
    !!
    !! Adds to a particle dungeon n particles, sampled
    !! from the corresponding source distributions
    !!
    !! Args:
    !!   dungeon [inout] -> particle dungeon to be added to
    !!   n [in]          -> number of particles to place in dungeon
    !!   rand [inout]    -> particle RNG object
    !!   matIdx [in]     -> optional unused argument, here so that subclasses can override to
    !!                      select matIdx to sample from
    !!
    !! Result:
    !!   A dungeon populated with n particles sampled from the source, plus
    !!   particles already present in dungeon
    !!
    subroutine append(self, dungeon, n, rand, matIdx)
      class(source), intent(inout)            :: self
      type(particleDungeon), intent(inout)    :: dungeon
      integer(shortInt), intent(in)           :: n
      class(RNG), intent(inout)               :: rand
      integer(shortInt), intent(in), optional :: matIdx
      integer(shortInt)                       :: i

      ! Generate n particles to populate dungeon
      do i = 1, n
        call dungeon % detain(self % sampleParticle(rand))
      end do

    end subroutine append

    !!
    !! Return to uninitialised state
    !!
    elemental subroutine kill(self)
      class(source), intent(inout) :: self

      self % geom => null()

    end subroutine kill

end module source_inter
