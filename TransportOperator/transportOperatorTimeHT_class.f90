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
    class(particle), intent(inout)               :: p
    type(tallyAdmin), intent(inout)              :: tally
    class(particleDungeon), intent(inout)        :: thisCycle
    class(particleDungeon), intent(inout)        :: nextCycle
    real(defReal)                                :: sigmaT, dTime, dColl
    logical(defBool)                             :: finished
    integer(shortInt)                            :: matIdx
    character(100), parameter :: Here = 'timeTracking (transportOperatorTimeHT_class.f90)' 

    finished = .false.

    ! Get majorant XS inverse: 1/Sigma_majorant
    self % majorant_inv = ONE / self % xsData % getMajorantXS(p)

    trackingLoop:do

      ! Check for errors
      if (p % time /= p % time) call fatalError(Here, 'Particle time is NaN')

      ! Obtain sigmaT
      sigmaT = self % xsData % getTransMatXS(p, p % matIdx())

      ! Find distance to time boundary
      dTime = p % getSpeed() * (p % timeMax - p % time)

      ! Sample distance to move particle before collision
      dColl = -log( p % pRNG % get() ) / sigmaT

      ! Decide whether to use delta tracking or surface tracking
      ! Vastly different opacities make delta tracking infeasable
      if(sigmaT * self % majorant_inv > ONE - self % cutoff) then
        ! Delta tracking
        call self % deltaTracking(p, dTime, dColl, finished)
      else
        ! Surface tracking
        call self % surfaceTracking(p, dTime, dColl, finished)
      end if

      ! Check for particle leakage
      if (p % matIdx() == OUTSIDE_FILL) then
        p % fate = LEAK_FATE
        p % isDead = .true.
        exit trackingLoop
      end if

      ! Exit if transport is finished
      if (finished .eqv. .true.) exit trackingLoop

    end do trackingLoop

    call tally % reportTrans(p)

  end subroutine timeTracking

  !!
  !! Perform surface tracking
  !!
  subroutine surfaceTracking(self, p, dTime, dColl, finished)
    class(transportOperatorTimeHT), intent(inout) :: self
    class(particle), intent(inout)                :: p
    real(defReal), intent(in)                     :: dTime
    real(defReal), intent(in)                     :: dColl
    logical(defBool), intent(inout)               :: finished
    real(defReal)                                 :: dist
    integer(shortInt)                             :: event
    character(100), parameter :: Here = 'surfaceTracking (transportOperatorTimeHT_class.f90)'

    dist = min(dTime, dColl)

    ! Move through geometry using minimum distance
    call self % geom % move(p % coords, dist, event)

    p % time = p % time + dist / p % getSpeed()

    ! Check result of transport
    if (dist == dTime) then
      ! Time boundary
      if (event /= COLL_EV) call fatalError(Here, 'Moving dTime should result in COLL_EV')
      p % fate = AGED_FATE
      if (abs(p % time - p % timeMax) > 0.000001) call fatalError(Here, 'Particle time incorrect?')
      p % time = p % timeMax
      finished = .true.
    else if (dist == dColl) then
      ! Collision, increase time accordingly
      if (event /= COLL_EV) call fatalError(Here, 'Moving dColl should result in COLL_EV')
      p % time = p % time + dColl / p % getSpeed()
      finished = .true.
    end if

  end subroutine surfaceTracking

  !!
  !! Perform delta tracking
  !!
  subroutine deltaTracking(self, p, dTime, dColl, finished)
    class(transportOperatorTimeHT), intent(inout) :: self
    class(particle), intent(inout)                :: p
    real(defReal), intent(in)                     :: dTime
    real(defReal), intent(in)                     :: dColl
    logical(defBool), intent(inout)               :: finished
    real(defReal)                                 :: sigmaT

    ! Determine which distance to move particle
    if (dColl < dTime) then
      ! Move partice to potential collision location
      call self % geom % teleport(p % coords, dColl)
      p % time = p % time + dColl / p % getSpeed()
    else
      ! Move particle to end of time step location
      call self % geom % teleport(p % coords, dTime)
      p % fate = AGED_FATE
      p % time = p % timeMax
      finished = .true.
      return
    end if

    ! Check for particle leakage
    if (p % matIdx() == OUTSIDE_FILL) return

    ! Obtain local cross-section
    sigmaT = self % xsData % getTransMatXS(p, p % matIdx())

    ! Roll RNG to determine if the collision is real or virtual
    ! Exit the loop if the collision is real
    if (p % pRNG % get() < sigmaT * self % majorant_inv) finished = .true.

  end subroutine deltaTracking

  !!
  !! Provide transport operator with delta tracking/surface tracking cutoff
  !!
  !! Cutoff of 1 gives exclusively delta tracking, cutoff of 0 gives exclusively surface tracking
  !!
  subroutine init(self, dict)
    class(transportOperatorTimeHT), intent(inout) :: self
    class(dictionary), intent(in)              :: dict
    class(dictionary), pointer                 :: tempDict
    integer(shortInt)                          :: nMats
    real(defReal), dimension(6)                :: bounds
    real(defReal)                              :: lengthScale

    ! Initialise superclass
    call init_super(self, dict)

    ! Get cutoff value
    call dict % getOrDefault(self % cutoff, 'cutoff', 0.7_defReal)

  end subroutine init

end module transportOperatorTimeHT_class
