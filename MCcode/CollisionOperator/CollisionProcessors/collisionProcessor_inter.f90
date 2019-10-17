module collisionProcessor_inter

  use numPrecision
  use endfConstants
  use genericProcedures,     only : fatalError, numToChar
  use dictionary_class,      only : dictionary
  use RNG_class,             only : RNG
  use particle_class,        only : particle
  use particleDungeon_class, only : particleDungeon

  ! Tally interfaces
  use tallyCodes
  use tallyAdmin_class,      only : tallyAdmin

  implicit none
  private

  !!
  !! Data package with all relevant data about the collision to move beetween customisable
  !! procedures
  !!
  type, public :: collisionData
    integer(shortInt) :: matIdx  = -1   !! Material Index at collision
    integer(shortInt) :: nucIdx  = -1   !! Nuclide Index of target
    integer(shortInt) :: MT      = 0    !! MT Number of Realction
    real(defReal)     :: muL     = ONE  !! Cosine of deflection angle in LAB-frame
    real(defReal)     :: A       = ZERO !! Target Mass [Neutron Mass]
    real(defReal)     :: kT      = ZERO !! Target temperature [MeV]
  end type


  !!
  !! This is an abstract interface for all types of collision processing
  !!  -> This interface only deals with SCALAR processing of collisions
  !!  -> Note that it is NOT just a collection of deferred function. There is a master
  !!     non_overridable function that controls the flow of the calculation and calls a number
  !!     of customisable deferred function.
  !!  -> Interesed user can refer to http://www.gotw.ca/publications/mill18.htm for justification
  !!     for this approach.
  !!
  !! Public interface:
  !!   collide(p, tally, thisCycle, nextCycle) -> Given particle, tallyAdmin, and particleDungeons
  !!     for particles in this cycle, or next cycle performs particle collision. Sends pre and post
  !!     collision events report to the tally admin. Reports end of history if particle was absorbed.
  !!
  !!   init(dict) -> initialises collisionProcessor from a dictionary
  !!
  !!
  !! Customisable procedures or collision actions (implemented in subclasses):
  !!   sampleCollision -> should determine collision target and type. Sets approperiate data in
  !!                      collisionData package.
  !!   implicit        -> preforms any implicit treatment to be done before reaction processing,
  !!                      e.g. implicit fission site production
  !!   scatter         -> defines behaviour for a scattering event (macroscopic scatter or any
  !!                      nuclide scattering reaction (elastic, inealastic, campton, NXN etc. )
  !!   capture         -> defines behaviour for any capture reaction (e.g. gamma capture, photoelectric)
  !!   fission         -> defines bahaviour for any fission reaction
  !!   cutoffs         -> Any post collision implicit treatments i.e. energy cutoffs
  !!
  type, public, abstract :: collisionProcessor
    private
  contains
    ! Master non-overridable procedures
    procedure, non_overridable :: collide

    ! Extendable initialisation procedure
    procedure :: init

    ! Customisable deffered procedures
    procedure(collisionAction),deferred  :: sampleCollision
    procedure(collisionAction),deferred  :: implicit
    procedure(collisionAction),deferred  :: elastic
    procedure(collisionAction),deferred  :: inelastic
    procedure(collisionAction),deferred  :: capture
    procedure(collisionAction),deferred  :: fission
    procedure(collisionAction),deferred  :: cutoffs

  end type collisionProcessor

  !! Extandable procedures
  public :: init


  abstract interface
    !!
    !! Procedure interface for all customisable actions associated with
    !! processing of sollision event (scatter, fission etc.)
    !!
    subroutine collisionAction(self, p, collDat, thisCycle, nextCycle)
      import :: collisionProcessor, &
                collisionData, &
                particle,&
                particleDungeon
      class(collisionProcessor), intent(inout) :: self
      class(particle), intent(inout)           :: p
      type(collisionData), intent(inout)       :: collDat
      class(particleDungeon),intent(inout)     :: thisCycle
      class(particleDungeon),intent(inout)     :: nextCycle
    end subroutine collisionAction
  end interface

contains

  !!
  !! Generic flow of collision processing
  !!
  subroutine collide(self, p, tally ,thisCycle, nextCycle)
    class(collisionProcessor), intent(inout) :: self
    class(particle), intent(inout)           :: p
    type(tallyAdmin), intent(inout)          :: tally
    class(particleDungeon),intent(inout)     :: thisCycle
    class(particleDungeon),intent(inout)     :: nextCycle
    type(collisionData)                      :: collDat
    character(100),parameter                 :: Here = ' collide (collisionProcessor.f90)'

    ! Load material index into data package
    collDat % matIdx = p % matIdx()

    ! Report in-collision & save pre-collison state
    call tally % reportInColl(p)
    call p % savePreCollision()

    ! Choose collision nuclide and general type (Scatter, Capture or Fission)
    call self % sampleCollision(p, collDat, thisCycle, nextCycle)

    ! Perform implicit treatment
    call self % implicit(p, collDat, thisCycle, nextCycle)

    ! Select physics to be processed based on MT number
    select case(collDat % MT)
      case(N_N_elastic, macroAllScatter)
        call self % elastic(p, collDat, thisCycle, nextCycle)

      case(N_N_inelastic)
        call self % inelastic(p, collDat, thisCycle, nextCycle)

      case(N_DISAP, macroCapture)
        call self % capture(p, collDat, thisCycle, nextCycle)

      case(N_FISSION, macroFission)
        call self % fission(p, collDat, thisCycle, nextCycle)

      case(noInteraction)
        ! Do nothing

      case default
        call fatalError(Here, 'Unsupported MT number: '// numToChar(collDat % MT))

    end select

    ! Apply post collision implicit treatments
    call self % cutoffs(p, collDat, thisCycle, nextCycle)

    ! Report out-of-collision
    call tally % reportOutColl(p, collDat % MT, collDat % muL)

    ! Report end-of-history if particle was killed
    if( p % isDead) then
      p % fate = ABS_FATE
      call tally % reportHist(p)
    end if

  end subroutine collide

  !!
  !! Extendable initialisation procedure
  !!
  subroutine init(self, dict)
    class(collisionProcessor), intent(inout) :: self
    class(dictionary), intent(in)            :: dict

    ! For now does nothing

  end subroutine init

end module collisionProcessor_inter
