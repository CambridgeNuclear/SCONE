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
  use geometry_inter,             only : geometry

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

  subroutine init(self, dict)
    class(transportOperatorHT), intent(inout) :: self
    class(dictionary), intent(in)             :: dict

    ! Initialise superclass
    call init_super(self, dict)

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
    majorant_inv = ONE / self % xsData % getTrackingXS(p, p % matIdx(), MAJORANT_XS)

    ! Obtain the local cross-section
    sigmaT = self % xsData % getTotalMatXS(p, p % matIdx())

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
    character(100), parameter :: Here = 'deltaTracking (transportOperatorHT_class.f90)'

    ! Get majorant XS inverse: 1/Sigma_majorant
    majorant_inv = ONE / self % xsData % getTrackingXS(p, p % matIdx(), MAJORANT_XS)

   ! Should never happen! Prevents Inf distances
    if (abs(majorant_inv) > huge(majorant_inv)) call fatalError(Here, "Majorant is 0")

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
      if(p % matIdx() == VOID_MAT) then
        call tally % reportInColl(p, .true.)
        cycle DTLoop
      end if

      ! Give error if the particle somehow ended in an undefined material
      if (p % matIdx() == UNDEF_MAT) then
        print *, p % rGlobal()
        call fatalError(Here, "Particle is in undefined material")
      end if

      ! Obtain the local cross-section
      sigmaT = self % xsData % getTotalMatXS(p, p % matIdx())

      ! Roll RNG to determine if the collision is real or virtual
      ! Exit the loop if the collision is real, report collision if virtual
      if (p % pRNG % get() < sigmaT*majorant_inv) then
        exit DTLoop
      else
        call tally % reportInColl(p, .true.)
      end if

    end do DTLoop

    call tally % reportTrans(p)
  end subroutine deltaTracking


  subroutine surfaceTracking(self, p, tally, thisCycle, nextCycle)
    class(transportOperatorHT), intent(inout) :: self
    class(particle), intent(inout)            :: p
    type(tallyAdmin), intent(inout)           :: tally
    class(particleDungeon),intent(inout)      :: thisCycle
    class(particleDungeon),intent(inout)      :: nextCycle
    integer(shortInt)                         :: event
    real(defReal)                             :: sigmaT, dist
    character(100), parameter :: Here = 'surfaceTracking (transportOperatorHT_class.f90)'

    STLoop: do

      ! Obtain the local cross-section
      if (p % matIdx() == VOID_MAT) then
        dist = INFINITY

      else
        sigmaT = self % xsData % getTrackingXS(p, p % matIdx(), MATERIAL_XS)
        dist = -log( p % pRNG % get()) / sigmaT

        ! Should never happen! Catches NaN distances
        if (dist /= dist) call fatalError(Here, "Distance is NaN")

      end if

      ! Save state before movement
      call p % savePrePath()

      ! Move to the next stop. NOTE: "move" resets dist to distanced moved!
      call self % geom % move(p % coords, dist, event)

      ! Send tally report for a path moved
      call tally % reportPath(p, dist)

      ! Kill particle if it has leaked
      if (p % matIdx() == OUTSIDE_FILL) then
        p % isDead = .true.
        p % fate = LEAK_FATE
      end if

      ! Give error if the particle somehow ended in an undefined material
      if (p % matIdx() == UNDEF_MAT) then
        print *, p % rGlobal()
        call fatalError(Here, "Particle is in undefined material")
      end if

      ! Return if particle stoped at collision (not cell boundary)
      if (event == COLL_EV .or. p % isDead) exit STLoop

    end do STLoop

    call tally % reportTrans(p)

  end subroutine surfaceTracking


end module transportOperatorHT_class
