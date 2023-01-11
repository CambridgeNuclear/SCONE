!!
!! Transport operator time-dependent problems using a hybrid of delta tracking and surface tracking
!!
module transportOperatorTimeHT_class
  use numPrecision
  use universalVariables

  use genericProcedures,          only : fatalError, numToChar
  use particle_class,             only : particle, P_PHOTON
  use particleDungeon_class,      only : particleDungeon
  use dictionary_class,           only : dictionary
  use rng_class,                  only : rng

  ! Superclass
  use transportOperator_inter,    only : transportOperator, init_super => init

  ! Tally interface
  use tallyCodes
  use tallyAdmin_class,           only : tallyAdmin

  ! Nuclear data interfaces
  use nuclearDatabase_inter,      only : nuclearDatabase

  ! Geometry interfaces
  use simpleGrid_class,           only : simpleGrid

  implicit none
  private

  !!
  !! Transport operator that moves a particle with using hybrid tracking, up to a time boundary
  !!
  type, public, extends(transportOperator)   :: transportOperatorTimeHT
    real(defReal)                            :: majorant_inv
    real(defReal)                            :: deltaT
    real(defReal)                            :: cutoff
  contains
    procedure          :: transit => timeTracking
    procedure          :: init
    procedure, private :: surfaceTracking
    procedure, private :: deltaTracking
  end type transportOperatorTimeHT

contains

  subroutine timeTracking(self, p, tally, thisCycle, nextCycle)
    class(transportOperatorTimeHT), intent(inout) :: self
    class(particle), intent(inout)                :: p
    type(tallyAdmin), intent(inout)               :: tally
    class(particleDungeon), intent(inout)         :: thisCycle
    class(particleDungeon), intent(inout)         :: nextCycle
    real(defReal)                                 :: sigmaT
    character(100), parameter :: Here = 'timeTracking (transportOperatorTimeHT_class.f90)' 

    ! Get majorant XS inverse: 1/Sigma_majorant
    if (associated(self % grid)) then
      self % majorant_inv = ONE / self % grid % getValue(p % coords % lvl(1) % r, p % coords % lvl(1) % dir)
    else
      self % majorant_inv = ONE / self % xsData % getMajorantXS(p)
    end if

    ! Check for errors
    if (p % time /= p % time) call fatalError(Here, 'Particle time is NaN')

    ! Obtain sigmaT
    sigmaT = self % xsData % getTransMatXS(p, p % matIdx())

    ! Decide whether to use delta tracking or surface tracking
    ! Vastly different opacities make delta tracking infeasable
    if(sigmaT * self % majorant_inv > ONE - self % cutoff) then
      ! Delta tracking
      call self % deltaTracking(p)
    else
      ! Surface tracking
      call self % surfaceTracking(p)
    end if

    ! Check for particle leakage
    if (p % matIdx() == OUTSIDE_FILL) then
      p % fate = LEAK_FATE
      p % isDead = .true.
    end if

    call tally % reportTrans(p)

  end subroutine timeTracking

  !!
  !! Perform surface tracking
  !!
  subroutine surfaceTracking(self, p)
    class(transportOperatorTimeHT), intent(inout) :: self
    class(particle), intent(inout)                :: p
    real(defReal)                                 :: dTime, dColl, dist, sigmaT
    integer(shortInt)                             :: event
    character(100), parameter :: Here = 'surfaceTracking (transportOperatorTimeHT_class.f90)'

    STLoop:do

      ! Find distance to time boundary
      dTime = lightSpeed * (p % timeMax - p % time)

      ! Sample distance to collision
      sigmaT = self % xsData % getTransMatXS(p, p % matIdx())
      dColl = -log( p % pRNG % get() ) / sigmaT

      ! Choose minimum distance
      dist = min(dTime, dColl)

      ! Move through geometry using minimum distance
      call self % geom % move(p % coords, dist, event)

      ! Check for particle leakage
      if (p % matIdx() == OUTSIDE_FILL) return

      ! Increase time based on distance moved
      p % time = p % time + dist / lightSpeed

      ! Check result of transport
      if (dist == dTime) then
        ! Time boundary
        if (event /= COLL_EV) call fatalError(Here, 'Move outcome should be COLL_EV &
                                                    &after moving dTime')
        p % fate = AGED_FATE
        p % time = p % timeMax
        exit STLoop

      else if (dist == dColl) then
        ! Collision, increase time accordingly
        if (event /= COLL_EV) call fatalError(Here, 'Move outcome should be COLL_EV &
                                                    &after moving dTime')
        exit STLoop

      end if

    end do STLoop

  end subroutine surfaceTracking

  !!
  !! Perform delta tracking
  !!
  subroutine deltaTracking(self, p)
    class(transportOperatorTimeHT), intent(inout) :: self
    class(particle), intent(inout)                :: p
    real(defReal)                                 :: dTime, dColl, dGrid, sigmaT
    character(100), parameter :: Here = 'deltaTracking (transportOperatorTimeHT_class.f90)'

    dColl = ZERO
    dGrid = INF
    if (associated(self % grid)) dGrid = self % grid % getDistance(p % coords % lvl(1) % r, p % coords % lvl(1) % dir)

    DTLoop:do

      dGrid = dGrid - dColl

      ! Find distance to time boundary
      dTime = lightSpeed * (p % timeMax - p % time)

      ! Sample distance to collision
      dColl = -log( p % pRNG % get() ) * self % majorant_inv

      if (dGrid < dTime .and. dGrid < dColl) then
        call self % geom % teleport(p % coords, dGrid)
        p % time = p % time + dGrid / lightSpeed
        self % majorant_inv = ONE / self % grid % getValue(p % coords % lvl(1) % r, p % coords % lvl(1) % dir)
        dGrid = self % grid % getDistance(p % coords % lvl(1) % r, p % coords % lvl(1) % dir)
        cycle DTLoop

      ! If dTime < dColl, move to end of time step location
      else if (dTime < dColl) then
        call self % geom % teleport(p % coords, dTime)
        p % fate = AGED_FATE
        p % time = p % timeMax
        exit DTLoop
      end if

      ! Otherwise, move to potential collision location
      call self % geom % teleport(p % coords, dColl)
      p % time = p % time + dColl / lightSpeed

      ! Check for particle leakage
      if (p % matIdx() == OUTSIDE_FILL) return

      ! Obtain local cross-section
      sigmaT = self % xsData % getTransMatXS(p, p % matIdx())

      ! Roll RNG to determine if the collision is real or virtual
      ! Exit the loop if the collision is real
      if (p % pRNG % get() < sigmaT * self % majorant_inv) exit DTLoop

      ! Protect against infinite loop
      if (sigmaT*self % majorant_inv == 0) call fatalError(Here, '100 % virtual collision chance, &
                                                                 &potentially infinite loop')

    end do DTLoop

  end subroutine deltaTracking

  !!
  !! Provide transport operator with delta tracking/surface tracking cutoff
  !!
  !! Cutoff of 1 gives exclusively delta tracking, cutoff of 0 gives exclusively surface tracking
  !!
  subroutine init(self, dict, grid)
    class(transportOperatorTimeHT), intent(inout)    :: self
    class(dictionary), intent(in)                    :: dict
    class(simpleGrid), intent(in), pointer, optional :: grid

    ! Initialise superclass
    call init_super(self, dict)

    ! Get cutoff value
    call dict % getOrDefault(self % cutoff, 'cutoff', 0.7_defReal)

    ! Store grid pointer
    if (present(grid)) self % grid => grid

  end subroutine init

end module transportOperatorTimeHT_class
