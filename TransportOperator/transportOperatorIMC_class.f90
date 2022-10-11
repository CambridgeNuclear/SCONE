!!
!! Transport operator for implicit Monte Carlo scheme using delta tracking
!!
module transportOperatorIMC_class
  use numPrecision
  use universalVariables

  use genericProcedures,          only : fatalError, numToChar
  use particle_class,             only : particle, P_PHOTON
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
  use mgIMCMaterial_inter,        only : mgIMCMaterial, mgIMCMaterial_CptrCast

  implicit none
  private

  !!
  !! Transport operator that moves a particle with IMC tracking
  !!
  type, public, extends(transportOperator) :: transportOperatorIMC
    class(mgIMCMaterial), pointer, public :: mat    => null()
    real(defReal)                         :: majorant_inv
  contains
    procedure          :: transit => imcTracking
    procedure, private :: materialTransform
    procedure, private :: surfaceTracking
    procedure, private :: deltaTracking
  end type transportOperatorIMC

contains

  subroutine imcTracking(self, p, tally, thisCycle, nextCycle)
    class(transportOperatorIMC), intent(inout) :: self
    class(particle), intent(inout)             :: p
    type(tallyAdmin), intent(inout)            :: tally
    class(particleDungeon), intent(inout)      :: thisCycle
    class(particleDungeon), intent(inout)      :: nextCycle
    real(defReal)                              :: sigmaT, dTime, dColl
    logical(defBool)                           :: finished
    character(100), parameter :: Here = 'IMCTracking (transportOperatorIMC_class.f90)' 

    finished = .false.

    ! Get majorant XS inverse: 1/Sigma_majorant
    self % majorant_inv = ONE / self % xsData % getMajorantXS(p)

    ! Deal with material particles, only relevant for ISMC
    if(p % getType() == P_MATERIAL_MG) then
      call self % materialTransform(p, tally)
      if(p % fate == TIME_FATE) return
    end if

    IMCLoop:do

      ! Check for errors
      if (p % getType() /= P_PHOTON_MG) call fatalError(Here, 'Particle is not MG Photon')
      if (p % time /= p % time) call fatalError(Here, 'Particle time is NaN')

      ! Obtain sigmaT
      sigmaT = self % xsData % getTransMatXS(p, p % matIdx())

      ! Find distance to time boundary
      dTime = lightSpeed * (p % timeMax - p % time)

      ! Sample distance to move particle before collision
      dColl = -log( p % pRNG % get() ) / sigmaT

      ! Decide whether to use delta tracking or surface tracking
      ! Vastly different opacities make delta tracking infeasable
      if(sigmaT * self % majorant_inv > 0.3) then
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
        exit IMCLoop
      end if

      ! Exit if transport is finished
      if (finished .eqv. .true.) exit IMCLoop

    end do IMCLoop

    call tally % reportTrans(p)

  end subroutine imcTracking


  !!
  !! Transform material particles into radiation photons with
  !! probability per unit time of c*sigma_a*fleck*eta
  !!
  !! Used only for ISMC, not for standard IMC
  !!
  subroutine materialTransform(self, p, tally)
    class(transportOperatorIMC), intent(inout) :: self
    class(particle), intent(inout)             :: p
    type(tallyAdmin), intent(inout)            :: tally
    real(defReal)                              :: sigmaT, fleck, eta, mu, phi
    real(defReal), dimension(3)                :: dir
    character(100), parameter                  :: Here = 'materialTransform (transportOperatorIMC_class.f90)'

    ! Confirm that time = 0
    !if (p % time .ne. 0) call fatalError(Here, 'Material particle should have time = 0')

    ! Get and verify material pointer
    self % mat => mgIMCMaterial_CptrCast( self % xsData % getMaterial( p % matIdx()))
    if(.not.associated(self % mat)) call fatalError(Here, "Failed to get MG IMC Material")

    sigmaT = self % xsData % getTransMatXS(p, p % matIdx())     !! Should be sigma_a, may need changing when sorting out cross-sections
    fleck  = self % mat % getFleck()
    eta    = self % mat % getEta()

    ! Sample time to transform into radiation photon
    p % time = p % time - log(p % pRNG % get()) / (sigmaT*fleck*eta*lightSpeed)

    ! Deal with eta = 0
    if (p % time /= p % time) p % time = INF

    ! Exit loop if particle remains material until end of time step
    if (p % time >= p % timeMax) then
      p % fate = TIME_FATE
      p % time = p % timeMax
      ! Tally energy for next temperature calculation
      call tally % reportHist(p)
    else
      p % type = P_PHOTON
      ! Resample direction
      mu = 2 * p % pRNG % get() - 1
      phi = p % pRNG % get() * 2*pi
      dir(1) = mu
      dir(2) = sqrt(1-mu**2) * cos(phi)
      dir(3) = sqrt(1-mu**2) * sin(phi)
      call p % point(dir)
    end if

  end subroutine materialTransform

  !!
  !! Perform surface tracking
  !!
  subroutine surfaceTracking(self, p, dTime, dColl, finished)
    class(transportOperatorIMC), intent(inout) :: self
    class(particle), intent(inout)             :: p
    real(defReal), intent(in)                  :: dTime
    real(defReal), intent(in)                  :: dColl
    logical(defBool), intent(inout)            :: finished
    real(defReal)                              :: dist
    integer(shortInt)                          :: event
    character(100), parameter                  :: Here = 'surfaceTracking (transportOperatorIMC_class.f90)'

    dist = min(dTime, dColl)

    ! Move through geometry using minimum distance
    call self % geom % move(p % coords, dist, event)

    p % time = p % time + dist / lightSpeed

    ! Check result of transport
    if (dist == dTime) then
      ! Time boundary
      if (event /= COLL_EV) call fatalError(Here, 'Move outcome should be COLL_EV after moving dTime')
      p % fate = TIME_FATE
      if (abs(p % time - p % timeMax)>0.000001) call fatalError(Here, 'Particle time is somehow incorrect')
      p % time = p % timeMax
      finished = .true.
    else if (dist == dColl) then
      ! Collision, increase time accordingly
      if (event /= COLL_EV) call fatalError(Here, 'Move outcome should be COLL_EV after moving dTime')
      finished = .true.
    end if

  end subroutine surfaceTracking

  !!
  !! Perform delta tracking
  !!
  subroutine deltaTracking(self, p, dTime, dColl, finished)
    class(transportOperatorIMC), intent(inout) :: self
    class(particle), intent(inout)             :: p
    real(defReal), intent(in)                  :: dTime
    real(defReal), intent(in)                  :: dColl
    logical(defBool), intent(inout)            :: finished
    real(defReal)                              :: sigmaT

    ! Determine which distance to move particle
    if (dColl < dTime) then
      ! Move partice to potential collision location
      call self % geom % teleport(p % coords, dColl)
      p % time = p % time + dColl / lightSpeed
    else
      ! Move particle to end of time step location
      call self % geom % teleport(p % coords, dTime)
      p % fate = TIME_FATE
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

end module transportOperatorIMC_class
