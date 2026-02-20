module transportOperator_inter

  use numPrecision
  use universalVariables
  use errors_mod,               only : fatalError

  use particle_class,           only : particle
  use particleDungeon_class,    only : particleDungeon
  use dictionary_class,         only : dictionary


  ! Geometry interfaces
  use geometryReg_mod,          only : gr_geomPtr => geomPtr, &
                                       gr_hasField => hasField, &
                                       gr_fieldIdx => fieldIdx, &
                                       gr_fieldPtr => fieldPtr
  use geometry_inter,           only : geometry
  
  use field_inter,              only : field
  use pieceConstantField_inter, only : pieceConstantField, pieceConstantField_CptrCast

  ! Tally interface
  use tallyAdmin_class,         only : tallyAdmin

  ! Nuclear data interfaces
  use nuclearDataReg_mod,       only : ndReg_get => get
  use nuclearDatabase_inter,    only : nuclearDatabase



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
  !! Procedures generic to transport operators:
  !!   localConditions(p) -> obtains local conditions of temperature and density for use in transport
  !!
  type, abstract, public :: transportOperator
    !! Nuclear Data block pointer -> public so it can be used by subclasses (protected member)
    class(nuclearDatabase), pointer :: xsData => null()

    !! Geometry pointer -> public so it can be used by subclasses (protected member)
    class(geometry), pointer         :: geom        => null()

  contains
    ! Public interface
    procedure, non_overridable :: transport

    ! Extentable initialisation and deconstruction procedure
    procedure :: init
    procedure :: kill

    ! Query for local conditions of temperature and density.
    procedure, non_overridable :: localConditions

    ! Customisable deferred procedures
    procedure(transit), deferred :: transit

  end type transportOperator

  ! Extendable procedures
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
    if (p  % isDead) then
      call tally % reportHist(p)
    end if

  end subroutine transport

  !!
  !! Queries local conditions of temperature and density from
  !! the geometry registry. Updates the particle to carry this info.
  !!
  subroutine localConditions(self, p)
    class(transportOperator), intent(in) :: self
    class(particle), intent(inout)       :: p
    class(field), pointer                :: genericField
    class(pieceConstantField), pointer   :: pcField
    
    ! Temperature check
    if (gr_hasField(nameTemperature)) then
      genericField => gr_fieldPtr(gr_fieldIdx(nameTemperature))
      pcField => pieceConstantField_CptrCast(genericField)
      p % T = pcField % at(p % coords)
    end if

    ! Density check
    if (gr_hasField(nameDensity)) then
      genericField => gr_fieldPtr(gr_fieldIdx(nameDensity))
      pcField => pieceConstantField_CptrCast(genericField)
      p % rho = pcField % at(p % coords)
    end if

  end subroutine localConditions

  !!
  !! Initialise transport operator from dictionary and geometry
  !!
  subroutine init(self, dict)
    class(transportOperator), intent(inout)  :: self
    class(dictionary), intent(in)            :: dict

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
