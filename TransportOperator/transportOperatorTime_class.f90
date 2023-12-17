!!
!! Surface tracking transport operator for time-dependent problems
!!
module transportOperatorTime_class
  use numPrecision
  use universalVariables

  use genericProcedures,          only : fatalError, numToChar
  use particle_class,             only : particle, P_PHOTON, P_MATERIAL
  use particleDungeon_class,      only : particleDungeon
  use dictionary_class,           only : dictionary
  use rng_class,                  only : rng

  ! Superclass
  use transportOperator_inter,    only : transportOperator

  ! Geometry interfaces
  use geometry_inter,             only : geometry

  ! Tally interface
  use tallyCodes
  use tallyAdmin_class,           only : tallyAdmin

  ! Nuclear data interfaces
  use nuclearDatabase_inter,      only : nuclearDatabase
  use mgIMCDatabase_inter,        only : mgIMCDatabase, mgIMCDatabase_CptrCast

  implicit none
  private

  !!
  !! Transport operator that moves a particle using surface tracking, up to a time boundary
  !!
  type, public, extends(transportOperator)   :: transportOperatorTime
  contains
    procedure          :: transit => timeTracking
    procedure, private :: surfaceTracking
    procedure, private :: materialTransform
  end type transportOperatorTime

contains

  subroutine timeTracking(self, p, tally, thisCycle, nextCycle)
    class(transportOperatorTime), intent(inout) :: self
    class(particle), intent(inout)              :: p
    type(tallyAdmin), intent(inout)             :: tally
    class(particleDungeon), intent(inout)       :: thisCycle
    class(particleDungeon), intent(inout)       :: nextCycle
    character(100), parameter :: Here = 'timeTracking (transportOperatorTime_class.f90)' 

    ! Transform material particles into photons
    if (p % type == P_MATERIAL) then
      call self % materialTransform(p, tally)
      ! Exit at time boundary
      if (p % fate == AGED_FATE) return
    end if

    ! Check for particle leakage
    if (p % matIdx() == OUTSIDE_FILL) then
      ! TODO: Figure out why this sometimes happens
      print *, 'WARNING: Leak before transport?'
      p % fate = LEAK_FATE
      p % isDead = .true.
      return
    end if

    call self % surfaceTracking(p)

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
    class(transportOperatorTime), intent(inout) :: self
    class(particle), intent(inout)              :: p
    real(defReal)                               :: dTime, dColl, dist, sigmaT
    integer(shortInt)                           :: event
    character(100), parameter :: Here = 'surfaceTracking (transportOperatorTime_class.f90)'

    STLoop:do

      ! Find distance to time boundary
      dTime = lightSpeed * (p % timeMax - p % time)

      ! Sample distance to collision
      sigmaT = self % xsData % getTransMatXS(p, p % matIdx())
      dColl = -log( p % pRNG % get() ) / sigmaT

      ! Ensure particle does not remain exactly on a boundary if dColl is close to 0
      if (event == CROSS_EV .and. dColl < SURF_TOL) then
        dColl = SURF_TOL
      end if

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
        p % fate = AGED_FATE
        p % time = p % timeMax
        exit STLoop

      else if (dist == dColl) then
        ! Collision
        exit STLoop

      end if

      if (event == COLL_EV) call fatalError(Here, 'Move outcome should be CROSS_EV or BOUNDARY_EV')

    end do STLoop

    if (event /= COLL_EV) call fatalError(Here, 'Move outcome should be COLL_EV')

  end subroutine surfaceTracking

  !!
  !! Determine when a material particle will transform into a photon for ISMC calculations
  !!
  !! Args:
  !!   p [inout]     -> material particle to be transformed
  !!   tally [inout] -> tally to keep track of material particles surviving time step
  !!
  subroutine materialTransform(self, p, tally)
    class(transportOperatorTime), intent(inout) :: self
    class(particle), intent(inout)              :: p
    type(tallyAdmin), intent(inout)             :: tally
    real(defReal)                               :: transformTime, mu, phi
    real(defReal), dimension(3)                 :: dir
    class(mgIMCDatabase), pointer               :: nucData
    integer(shortInt)                           :: matIdx, uniqueID
    character(100), parameter :: Here = 'materialTransform (transportOperatorIMC_class.f90)'

    ! Get pointer to nuclear database
    nucData => mgIMCDatabase_CptrCast(self % xsData)
    if (.not. associated(nucData)) call fatalError(Here, 'Unable to find mgIMCDatabase')

    ! Material particles can occasionally have coords placed in void if within SURF_TOL of boundary
    matIdx = p % matIdx()
    ! If so, get matIdx based on exact position (no adjustment for surface tol)
    if (matIdx == VOID_MAT .or. matIdx == OUTSIDE_MAT) then
      call self % geom % whatIsAt(matIdx, uniqueID, p % coords % lvl(1) % r, [ZERO,ZERO,ZERO])
    end if
    ! If still in invalid region, call fatalError
    if (matIdx == 0)        call fatalError(Here, 'Outside material particle')
    if (matIdx == VOID_MAT) call fatalError(Here, 'Void material particle')

    ! Sample time until emission
    transformTime = nucData % sampleTransformTime(matIdx, p % pRNG)
    p % time = min(p % timeMax, p % time + transformTime)

    ! Exit loop if particle remains material until end of time step
    if (p % time == p % timeMax) then
      p % fate = AGED_FATE
      ! Tally energy for next temperature calculation
      call tally % reportHist(p)

    ! Transform into photon
    else
      p % type = P_PHOTON
      ! Resample direction
      mu = 2 * p % pRNG % get() - 1
      phi = p % pRNG % get() * 2*pi
      dir(1) = mu
      dir(2) = sqrt(1-mu**2) * cos(phi)
      dir(3) = sqrt(1-mu**2) * sin(phi)
      call p % point(dir)
      ! Resample energy
      p % G = nucData % sampleEnergyGroup(matIdx, p % pRNG)
    end if

  end subroutine materialTransform

end module transportOperatorTime_class
