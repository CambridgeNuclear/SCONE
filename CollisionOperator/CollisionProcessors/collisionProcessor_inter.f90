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
    real(defReal)     :: E       = ZERO !! Collision energy (could be relative to target) [MeV]
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
  !!   alphaProd       -> defines behaviour for time production reactions
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
    procedure(collisionAction),deferred  :: alphaProd

  end type collisionProcessor

  !! Extandable procedures
  public :: init


  abstract interface
    !!
    !! Procedure interface for all customisable actions associated with
    !! processing of sollision event (scatter, fission etc.)
    !!
    subroutine collisionAction(self, p, tally, collDat, thisCycle, nextCycle)
      import :: collisionProcessor, &
                collisionData, &
                tallyAdmin, &
                particle,&
                particleDungeon
      class(collisionProcessor), intent(inout) :: self
      class(particle), intent(inout)           :: p
      type(tallyAdmin), intent(inout)          :: tally
      type(collisionData), intent(inout)       :: collDat
      class(particleDungeon),intent(inout)     :: thisCycle
      class(particleDungeon),intent(inout)     :: nextCycle
    end subroutine collisionAction
  end interface

contains

  !!
  !! Generic flow of collision processing
  !!
  subroutine collide(self, p, tally, thisCycle, nextCycle)
    class(collisionProcessor), intent(inout) :: self
    class(particle), intent(inout)           :: p
    type(tallyAdmin), intent(inout)          :: tally
    class(particleDungeon),intent(inout)     :: thisCycle
    class(particleDungeon),intent(inout)     :: nextCycle
    type(collisionData)                      :: collDat
    logical(defBool)                         :: virtual
    integer(shortInt)                        :: addCollision
    character(100),parameter                 :: Here = 'collide (collisionProcessor.f90)'

    ! Load material index into data package
    collDat % matIdx = p % matIdx()

    ! Choose collision nuclide and general type (Scatter, Capture or Fission)
    call self % sampleCollision(p, tally, collDat, thisCycle, nextCycle)
    
    ! In case of a TMS rejection, set collision as virtual
    if (collDat % MT == noInteraction) then
      virtual = .true.
      addCollision = 0
    else
      virtual = .false.
      addCollision = 1
    end if

    ! Report in-collision & save pre-collison state
    ! Note: the ordering must not be changed between feeding the particle to the tally
    ! and updating the particle's preCollision state, otherwise this may cause certain
    ! tallies (e.g., collisionProbability) to return dubious results
    call tally % reportInColl(p, virtual)

    call p % savePreCollision()

    ! Perform implicit treatment
    ! Don't include interactions which avoid sampling a nuclide
    if (collDat % MT /= noInteraction .and. collDat % MT /= N_TIME_ABS &
            .and. collDat % MT /= N_TIME_PROD) then
      call self % implicit(p, tally, collDat, thisCycle, nextCycle)
    end if

    ! Select physics to be processed based on MT number
    select case(collDat % MT)
      case(N_N_elastic, macroAllScatter)
        call self % elastic(p, tally, collDat, thisCycle, nextCycle)

      case(N_N_inelastic, macroIEScatter)
        call self % inelastic(p, tally, collDat, thisCycle, nextCycle)

      case(N_DISAP, macroDisappearance, N_TIME_ABS)
        call self % capture(p, tally, collDat, thisCycle, nextCycle)

      case(N_FISSION, macroFission)
        call self % fission(p, tally, collDat, thisCycle, nextCycle)

      case(N_TIME_PROD)
        call self % alphaProd(p, tally, collDat, thisCycle, nextCycle)

      case(noInteraction)
        ! Do nothing

      case default
        call fatalError(Here, 'Unsupported MT number: '// numToChar(collDat % MT))

    end select

    ! Apply post collision implicit treatments
    call self % cutoffs(p, tally, collDat, thisCycle, nextCycle)

    ! Update particle collision counter
    p % collisionN = p % collisionN + addCollision

    ! Report out-of-collision
    call tally % reportOutColl(p, collDat % MT, collDat % muL)

    ! Report end-of-history if particle was killed
    if (p % isDead) then
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
