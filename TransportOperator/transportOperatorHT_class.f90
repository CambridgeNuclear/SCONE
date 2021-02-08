!!
!! Transport operator for hybrid tracking
!!
module transportOperatorHT_class
  use numPrecision
  use universalVariables

  use genericProcedures,          only : fatalError, numToChar
  use particle_class,             only : particle
  use particleDungeon_class,      only : particleDungeon
  use dictionary_class,           only : dictionary
  use rng_class,                  only : rng

  ! Tally interface
  use tallyCodes
  use tallyAdmin_class,           only : tallyAdmin

  ! Superclass
  use transportOperator_inter,    only : transportOperator, init_super => init

  ! Geometry interfaces
  use cellGeometry_inter,         only : cellGeometry

  ! Nuclear data interfaces
  use nuclearDatabase_inter,      only : nuclearDatabase

  implicit none
  private

  !!
  !! Transport operator that moves a particle with hybrid tracking
  !!
  type, public, extends(transportOperator) :: transportOperatorHT
    !! Cutoff threshold between ST and DT
    real(defReal)             :: cutoff
  contains
    procedure :: init
    procedure :: transit => tracking_selection
    procedure, private :: deltaTracking
    procedure, private :: surfaceTracking
  end type transportOperatorHT

contains

  subroutine init(self, dict, geom)
    class(transportOperatorHT), intent(inout)  :: self
    class(dictionary), intent(in)              :: dict
    class(cellGeometry), pointer, intent(in)   :: geom

    ! Initialise superclass
    call init_super(self, dict, geom)

    ! Initialise this class
    call dict % getOrDefault(self % cutoff,'cutoff',0.9_defReal)

  end subroutine init

  subroutine tracking_selection(self, p, tally, thisCycle, nextCycle)
    class(transportOperatorHT), intent(inout)              :: self
    class(particle), intent(inout)                         :: p
    type(tallyAdmin), intent(inout)                        :: tally
    class(particleDungeon), intent(inout)                  :: thisCycle
    class(particleDungeon), intent(inout)                  :: nextCycle
    real(defReal)                                          :: majorant_inv, sigmaT, ratio
    character(100), parameter :: Here = 'hybridTracking (transportOIperatorHT_class.f90)'

    ! Get majornat XS inverse: 1/Sigma_majorant
    majorant_inv = ONE / self % xsData % getMajorantXS(p)

    ! Obtain the local cross-section
    sigmaT = self % xsData % getTransMatXS(p, p % matIdx())

    ! Calculate ratio between local cross-section and majorant
    ratio = sigmaT*majorant_inv

    ! Cut-off criterion to decide on tracking method
    if (ratio > (ONE - self % cutoff)) then
      call deltaTracking(self, p, tally, thisCycle, nextCycle)
    else
      call surfaceTracking(self, p, tally, thisCycle, nextCycle)
    end if

  end subroutine tracking_selection


  subroutine deltaTracking(self, p, tally, thisCycle, nextCycle)
    class(transportOperatorHT), intent(inout) :: self
    class(particle), intent(inout)            :: p
    type(tallyAdmin), intent(inout)           :: tally
    class(particleDungeon), intent(inout)     :: thisCycle
    class(particleDungeon), intent(inout)     :: nextCycle
    real(defReal)                             :: majorant_inv, sigmaT, distance
    character(100), parameter :: Here = 'deltaTracking (transportOperatorDT_class.f90)'

    ! Get majornat XS inverse: 1/Sigma_majorant
    majorant_inv = ONE / self % xsData % getMajorantXS(p)

    DTLoop:do
      distance = -log( p% pRNG % get() ) * majorant_inv

      ! Move partice in the geometry
      call self % geom % teleport(p % coords, distance)

      ! If particle has leaked exit
      if (p % matIdx() == OUTSIDE_FILL) then
        p % fate = LEAK_FATE
        p % isDead = .true.
        return
      end if

      ! Check for void
      if( p % matIdx() == VOID_MAT) cycle DTLoop

      ! Obtain the local cross-section
      sigmaT = self % xsData % getTransMatXS(p, p % matIdx())

      ! Roll RNG to determine if the collision is real or virtual
      ! Exit the loop if the collision is real
      if (p % pRNG % get() < sigmaT*majorant_inv) exit DTLoop

    end do DTLoop

    call tally % reportTrans(p)
  end subroutine deltaTracking


  subroutine surfaceTracking(self, p, tally, thisCycle, nextCycle)
    class(transportOperatorHT), intent(inout) :: self
    class(particle), intent(inout)            :: p
    type(tallyAdmin), intent(inout)           :: tally
    class(particleDungeon),intent(inout)      :: thisCycle
    class(particleDungeon),intent(inout)      :: nextCycle
    logical(defBool)                          :: isColl
    real(defReal)                             :: sigmaT, dist

    STLoop: do

      ! Obtain the local cross-section
      if( p % matIdx() == VOID_MAT) then
        dist = INFINITY

      else
        sigmaT = self % xsData % getTransMatXS(p, p % matIdx())
        dist = -log( p % pRNG % get()) / sigmaT
      end if

      ! Save state before movement
      call p % savePrePath()

      ! Move to the next stop. NOTE: "move" resets dist to distanced moved!
      call self % geom % move(p % coords, dist, isColl)

      ! Send tally report for a path moved
      call tally % reportPath(p, dist)

      ! Kill particle if it has leaked
      if( p % matIdx() == OUTSIDE_FILL) then
        p % isDead = .true.
        p % fate = LEAK_FATE
      end if

      ! Return if particle stoped at collision (not cell boundary)
      if( isColl .or. p % isDead) exit STLoop

    end do STLoop

    call tally % reportTrans(p)

  end subroutine surfaceTracking


end module transportOperatorHT_class
