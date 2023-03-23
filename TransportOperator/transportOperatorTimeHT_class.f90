!!
!! Transport operator for time-dependent problems using a hybrid of delta tracking and surface tracking
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
  !! Tracking method
  !!
  integer(shortInt), parameter :: HT = 1    ! Hybrid tracking
  integer(shortInt), parameter :: GT = 2    ! Grid tracking
  integer(shortInt), parameter :: ST = 3    ! Surface tracking
  integer(shortInt), parameter :: DT = 4    ! Delta tracking

  !!
  !! Transport operator that moves a particle with using hybrid tracking, up to a time boundary
  !!
  type, public, extends(transportOperator)   :: transportOperatorTimeHT
    real(defReal)                            :: deltaT
    real(defReal)                            :: cutoff
    integer(shortInt)                        :: method
  contains
    procedure          :: transit => timeTracking
    procedure          :: init
    procedure, private :: surfaceTracking
    procedure, private :: deltaTracking
    procedure, private :: getMajInv
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

    ! Select action based on specified method - HT and GT start with DT but can switch to ST
    if (self % method == ST) then
      call self % surfaceTracking(p)
    else
      call self % deltaTracking(p)
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


      ! TODO: Option to switch back to DT?


    end do STLoop

  end subroutine surfaceTracking

  !!
  !! Perform delta tracking - option to switch to surface tracking for HT and GT methods
  !!
  subroutine deltaTracking(self, p)
    class(transportOperatorTimeHT), intent(inout) :: self
    class(particle), intent(inout)                :: p
    real(defReal)                                 :: dTime, dColl, dGrid, sigmaT, majorant_inv, dist
    character(100), parameter :: Here = 'deltaTracking (transportOperatorTimeHT_class.f90)'

    ! Get majorant and grid crossing distance if required
    majorant_inv = self % getMajInv(p)
    if (self % method == GT) then
      dGrid = self % grid % getDistance(p % coords % lvl(1) % r, p % coords % lvl(1) % dir)
    else
      dGrid = INF
    end if

    ! Get initial opacity
    sigmaT = self % xsData % getTransMatXS(p, p % matIdx())

    DTLoop:do

      ! Switch to ST if required
      if (self % method /= DT) then
        if (sigmaT * majorant_inv < ONE - self % cutoff) then
          call self % surfaceTracking(p)
          return
        end if
      end if

      ! Find distance to time boundary
      dTime = lightSpeed * (p % timeMax - p % time)

      ! Sample distance to collision
      dColl = -log( p % pRNG % get() ) * majorant_inv

      ! Select particle by minimum distance
      dist = min(dColl, dTime, dGrid)
      call self % geom % teleport(p % coords, dist)
      p % time = p % time + dist / lightSpeed

      ! Exit in the case of particle leakage
      if (p % matIdx() == OUTSIDE_FILL) then
        return
      end if

      ! Act based on distance moved
      if (dist == dGrid) then
        ! Update values and cycle loop
        majorant_inv = self % getMajInv(p)
        dGrid = self % grid % getDistance(p % coords % lvl(1) % r, p % coords % lvl(1) % dir)
        sigmaT = self % xsData % getTransMatXS(p, p % matIdx())
        cycle DTLoop

      else if (dist == dTime) then
        ! Update particle fate and exit
        p % fate = AGED_FATE
        p % time = p % timeMax
        exit DTLoop

      else ! Dist == dColl
        ! Check for real or virtual collision
        sigmaT = self % xsData % getTransMatXS(p, p % matIdx())
        if (p % pRNG % get() < sigmaT * majorant_inv) exit DTLoop
        ! Update grid distance
        dGrid = dGrid - dColl

      end if

    end do DTLoop

  end subroutine deltaTracking

  !!
  !!
  !!
  function getMajInv(self, p) result (maj_inv)
    class(transportOperatorTimeHT), intent(in) :: self
    class(particle), intent(in)                :: p
    real(defReal)                              :: maj_inv

    if (self % method == GT) then
      maj_inv = ONE / self % grid % getValue(p % coords % lvl(1) % r, p % coords % lvl(1) % dir)
    else
      maj_inv = ONE / self % xsData % getMajorantXS(p)
    end if

  end function getMajInv


  !!
  !! Provide transport operator with delta tracking/surface tracking cutoff
  !!
  !! Cutoff of 1 gives exclusively delta tracking, cutoff of 0 gives exclusively surface tracking
  !!
  subroutine init(self, dict)
    class(transportOperatorTimeHT), intent(inout)    :: self
    class(dictionary), intent(in)                    :: dict
!    class(simpleGrid), intent(in), pointer, optional :: grid
    character(nameLen)                               :: method
    class(dictionary),pointer                        :: tempdict
    character(100), parameter :: Here = "init (transportOperatorTimeHT_class.f90)"

    ! Initialise superclass
    call init_super(self, dict)

    ! Get tracking method
    call dict % getOrDefault(method, 'method', 'HT')

    select case (method)

      ! Hybrid tracking
      case ('HT')
        self % method = HT
        ! Get cutoff value
        call dict % get(self % cutoff, 'cutoff')

      ! Grid tracking
      case ('GT')
        self % method = GT
        ! Get cutoff value
        call dict % get(self % cutoff, 'cutoff')

        ! Initialise grid for hybrid tracking
        tempDict => dict % getDictPtr('grid')
        allocate(self % grid)
        call self % grid % init(tempDict)

      ! Surface tracking
      case ('ST')
        self % method = ST

      ! Delta tracking
      case ('DT')
        self % method = DT

      case default
        call fatalError(Here, 'Invalid tracking method given. Must be HT, ST or DT.')

    end select

  end subroutine init

end module transportOperatorTimeHT_class
