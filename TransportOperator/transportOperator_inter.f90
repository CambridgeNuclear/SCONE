module transportOperator_inter

  use numPrecision
  use universalVariables
  use genericProcedures,          only : fatalError

  use particle_class,             only : particle
  use particleDungeon_class,      only : particleDungeon
  use dictionary_class,           only : dictionary
  use RNG_class,                  only : RNG

  ! Geometry interfaces
  use geometryReg_mod,            only : gr_geomPtr => geomPtr
  use geometry_inter,             only : geometry

  ! Tally interface
  use tallyAdmin_class,           only : tallyAdmin

  ! Nuclear data interfaces
  use nuclearDataReg_mod,         only : ndReg_get => get
  use nuclearDatabase_inter,      only : nuclearDatabase

  ! Geometry interfaces
  use trackingGrid_class,         only : trackingGrid


  implicit none
  private


  !!
  !! This is an abstract interface for all types of transport processing
  !!   -> This interface only deals with scalar processing of particle transport
  !!   -> Assumes that particle moves without any external forces (assumes that particle
  !!      moves along straight lines between colisions)
  !!
  !! Public interface:
  !!   transport(p, tally, thisCycle, nextCycle) -> given particle, tally and particle dungeons
  !!     for particles in this and next cycle performs movement of a particle in the geometry.
  !!     Sends transistion report to the tally. Sends history report as well if particle dies.
  !!   init(dict, geom) -> initialises transport operator from a dictionary and pointer to a
  !!                       geometry
  !!
  !! Customisable procedures or transport actions
  !!   transit(p, tally, thisCycle, nextCycle) -> implements movement from collision to collision
  !!
  type, abstract, public :: transportOperator
    !! Nuclear Data block pointer -> public so it can be used by subclasses (protected member)
    class(nuclearDatabase), pointer :: xsData => null()

    !! Geometry pointer -> public so it can be used by subclasses (protected member)
    class(geometry), pointer         :: geom        => null()

    class(trackingGrid), pointer     :: grid => null()

  contains
    ! Public interface
    procedure, non_overridable :: transport

    ! Extentable initialisation and deconstruction procedure
    procedure :: init
    procedure :: kill

    ! Customisable deferred procedures
    procedure(transit), deferred :: transit

  end type transportOperator

  ! Extandable procedures
  public :: init
  public :: kill


  abstract interface
    !!
    !! Move particle from collision to collision
    !!  Kill particle if needed
    !!
    subroutine transit(self, p, tally, thisCycle, nextCycle)
      import :: transportOperator, &
                particle, &
                tallyAdmin, &
                particleDungeon
      class(transportOperator), intent(inout) :: self
      class(particle), intent(inout)          :: p
      type(tallyAdmin), intent(inout)         :: tally
      class(particleDungeon), intent(inout)   :: thisCycle
      class(particleDungeon), intent(inout)   :: nextCycle
    end subroutine transit
  end interface

contains

  !!
  !! Master non-overridable subroutine to perform transport
  !!  Performs everything common to all types of transport
  !!
  subroutine transport(self, p, tally, thisCycle, nextCycle)
    class(transportOperator), intent(inout) :: self
    class(particle), intent(inout)          :: p
    type(tallyAdmin), intent(inout)         :: tally
    class(particleDungeon), intent(inout)   :: thisCycle
    class(particleDungeon), intent(inout)   :: nextCycle
    character(100),parameter :: Here ='transport (transportOperator_inter.f90)'

    ! Get nuclear data pointer form the particle
    self % xsData => ndReg_get(p % getType())

    ! Save geometry pointer
    self % geom => gr_geomPtr(p % geomIdx)

    ! Save pre-transition state
    call p % savePreTransition()

    ! Perform transit
    call self % transit(p, tally, thisCycle, nextCycle)

    ! Send history reports if particle died
    if( p  % isDead) then
      call tally % reportHist(p)
    end if

  end subroutine transport

  !!
  !! Initialise transport operator from dictionary and geometry
  !!
  subroutine init(self, dict)
    class(transportOperator), intent(inout)    :: self
    class(dictionary), intent(in)              :: dict

    ! Do nothing

  end subroutine init

  !!
  !! Free memory. Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(transportOperator), intent(inout) :: self

    self % geom   => null()
    self % xsData => null()

  end subroutine kill

end module transportOperator_inter
