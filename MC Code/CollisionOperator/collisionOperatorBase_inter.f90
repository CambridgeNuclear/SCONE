module collisionOperatorBase_inter

  use numPrecision
  use endfConstants
  use genericProcedures,     only : fatalError
  use dictionary_class,      only : dictionary
  use RNG_class,             only : RNG
  use particle_class,        only : particle, phaseCoord
  use particleDungeon_class, only : particleDungeon

  ! Nuclear Data interface
  use nuclearData_inter,     only : nuclearData

  ! Tally interfaces
  use tallyAdminBase_class,  only : tallyAdminBase

  implicit none
  private

  !!
  !! Abstract interface of collision operator
  !! Provides variables that are likly to be shared by all collision operators
  !! Provides generic flow of collision oprocessing so derived collision operators only need
  !! to provide implementation of a specific physics processing i.e. of a capture or scattering
  !!
  !! Assumptions:
  !! sampleCollision -> sets self % nucIdx if data is per nuclide
  !!                 -> sets self % MT to one of the lumped reaction channels from endfConstants:
  !!                    anyScatter, anyCapture or anyFission
  !! scatter  -> needs to choose whether to call elastic, inelastic, or N_XN
  !!
  type, public,abstract :: collisionOperatorBase
    ! Part of the interface
    class(tallyAdminBase), pointer :: tally

    !* Colision Operator Local variables. Used by all collision operators
    !* Should be protected (C++ terminology), but it cannot be implemented with Fortran.
    !* Components are public but are not part of the interface (DO NOT ACCESS THEM DIRECTLY)
    class(RNG), pointer            :: pRNG   => null()  !! Pointer to random number generator
    integer(shortInt)              :: nucIdx = 0        !! Collision nuclide index
    integer(shortInt)              :: matIdx = 0        !! Current material index
    integer(shortInt)              :: MT     = 0        !! Collision reaction channel
    real(defReal)                  :: muL    = 0.0      !! Deflection angle in current collision
    integer(longInt)               :: collCount = 0



  contains
    ! Following procedures specify generic flow of collision processing
    procedure     :: collide

    ! Following procedures are defined in a specific implementation.
    procedure(collisionMacro),deferred  :: sampleCollision  !! Samples lumped reaction channel i.e. anyScatter
    procedure(collisionMacro),deferred  :: implicit         !! Performs implicit treatment of reactions i.e. generate fission sites
    procedure(collisionMacro),deferred  :: scatter          !! Translates between anyScatter and elastic, inelastic, N_XN
!    procedure(collisionMacro),deferred  :: elastic          !! Elastic collision processing
!    procedure(collisionMacro),deferred  :: inelastic        !! Inelastic collision processing
!    procedure(collisionMacro),deferred  :: N_XN             !! (n,Xn) colission processing
    procedure(collisionMacro),deferred  :: capture          !! Capture reaction processing
    procedure(collisionMacro),deferred  :: fission          !! Fission reaction processing
    procedure(collisionMacro),deferred  :: cutoffs          !! Apply energy and weight cutoffs
    procedure(init)          ,deferred  :: init             !! Initialise
  end type collisionOperatorBase

  abstract interface

    !!
    !! Abstract interface for action of a collision operator on a particle
    !! Many subroutines share the same interface:
    !!
    subroutine collisionMacro(self,p,thisCycle,nextCycle)
      import :: collisionOperatorBase, &
                particle,&
                particleDungeon
      class(collisionOperatorBase), intent(inout) :: self
      class(particle), intent(inout)              :: p
      class(particleDungeon),intent(inout)        :: thisCycle
      class(particleDungeon),intent(inout)        :: nextCycle

    end subroutine collisionMacro

    !!
    !! Interface of initialisation subroutine
    !!
    subroutine init(self,nucData,settings)
      import :: collisionOperatorBase, &
                nuclearData, &
                dictionary
      class(collisionOperatorBase), intent(inout) :: self
      class(nuclearData), pointer, intent(in)     :: nucData
      class(dictionary), intent(in)               :: settings
    end subroutine init
  end interface

  contains

    !!
    !! Generic flow of collision processing
    !! Does not need to be repeted in specific implementations of the interface
    !! Can be overloaded if it is required or for performance (explicit inlining)
    !!
    subroutine collide(self,p,thisCycle,nextCycle)
      class(collisionOperatorBase), intent(inout) :: self
      class(particle), intent(inout)              :: p
      class(particleDungeon),intent(inout)        :: thisCycle
      class(particleDungeon),intent(inout)        :: nextCycle
      type(phaseCoord)                            :: preState
      character(100),parameter                    :: Here = ' collide (collisionOperatorBase.f90)'

      ! Load particle data
      self % collCount = self % collCount + 1
      self % matIdx = p % matIdx()
      self % pRNG => p % pRNG
      self % nucIdx = 0
      self % MT = 0

      ! Report in-collision & save pre-collison state
      call self % tally % reportInColl(p)
      preState = p

      call self % sampleCollision(p,thisCycle,nextCycle)
      call self % implicit(p,thisCycle,nextCycle)

      ! Select physics to be processed based on MT number
      select case(self % MT)
        case(anyScatter, macroAllScatter)
          call self % scatter(p,thisCycle,nextCycle)

        case(anyCapture, macroCapture)
          call self % capture(p,thisCycle,nextCycle)

        case(anyFission, macroFission)
          call self % fission(p,thisCycle,nextCycle)

        case default
          call fatalError(Here, 'Unsupported MT number')

      end select

      call self % cutoffs(p,thisCycle,nextCycle)

      ! Report out-of-collision
      call self % tally % reportOutColl(preState, p, self%MT, self % mul)
    end subroutine collide

end module collisionOperatorBase_inter
