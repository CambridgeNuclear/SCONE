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
  use transportOperator_inter,    only : transportOperator, init_super => init

  ! Geometry interfaces
  use geometry_inter,             only : geometry

  ! Tally interface
  use tallyCodes
  use tallyAdmin_class,           only : tallyAdmin

  ! Nuclear data interfaces
  use nuclearDatabase_inter,      only : nuclearDatabase
  use mgIMCMaterial_inter,        only : mgIMCMaterial, mgIMCMaterial_CptrCast
  use materialMenu_mod,           only : mm_nMat => nMat

  implicit none
  private

  !!
  !! Transport operator that moves a particle with IMC tracking
  !!
  type, public, extends(transportOperator)   :: transportOperatorIMC
    class(mgIMCMaterial), pointer, public    :: mat    => null()
    real(defReal)                            :: majorant_inv
    real(defReal)                            :: deltaT
    real(defReal)                            :: cutoff
    real(defReal), dimension(:), allocatable :: matMajs
    real(defReal), dimension(:), allocatable :: sigmaLocal
    integer(shortInt), dimension(:,:), allocatable :: matConnections
    integer(shortInt)                        :: majMapN = 0
    real(defReal), dimension(3)              :: top     = ZERO
    real(defReal), dimension(3)              :: bottom  = ZERO
    integer(shortInt)                        :: steps
  contains
    procedure          :: transit => imcTracking
    procedure          :: init
    procedure          :: buildMajMap
    procedure          :: updateMajorants
    procedure, private :: materialTransform
    procedure, private :: surfaceTracking
    procedure, private :: deltaTracking
    procedure, private :: simpleParticle
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

    ! Deal with material particles, only relevant for ISMC
    if(p % getType() == P_MATERIAL_MG) then
      call self % materialTransform(p, tally)
      if(p % fate == TIME_FATE) return
    end if

    ! Get majorant for particle
    if (allocated(self % matMajs)) then
      self % majorant_inv = ONE / self % matMajs(p % matIdx())
    else
      self % majorant_inv = ONE / self % xsData % getMajorantXS(p)
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
      if(sigmaT * self % majorant_inv > self % cutoff) then
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
  !! Generate particles and follow them to attempt to reduce majorant opacity seen by each cell.
  !!
  !! Particles sampled uniformly within geometry, with random directions, then incrementally moved
  !! up to distance dTime. Moving increments are determined by the number of steps to be taken,
  !! either given directly in input or calculated from dTime and a given cell lengthscale.
  !!
  !! Majorant of each cell is set to be the maximum opacity seen by any particle originating in
  !! that cell.
  !!
  subroutine buildMajMap(self, rand, xsData)
    class(transportOperatorIMC), intent(inout) :: self
    class(RNG), intent(inout)                  :: rand
    class(nuclearDatabase), intent(in), pointer :: xsData
    type(particle)                             :: p
    integer(shortInt)                          :: i, j, matIdx
    real(defReal)                              :: mu, phi, dist
    logical(defBool)                           :: finished = .false.

    ! Check that subroutine should be called
    if (.not. allocated(self % matMajs)) return

    ! Point to nuclear data, as otherwise cannot access until after first particle transport
    self % xsData => xsData

    ! Reset array
    self % matConnections = 0

    ! Calculate distance increments
    dist = self % deltaT * lightSpeed / self % steps

    do i = 1, self % majMapN

      ! Sample particle
      call self % simpleParticle(p, rand)
      matIdx = p % matIdx()

      ! Incrementally transport particle up to a distance dTime
      do j = 1, self % steps

        call self % geom % teleport(p % coords, dist)
        if (p % matIdx() == VOID_MAT .or. p % matIdx() == OUTSIDE_MAT) exit

        ! Update matConnections to signify a connection between starting mat and new mat
        self % matConnections(matIdx, p % matIdx()) = 1
        self % matConnections(p % matIdx(), matIdx) = 1

      end do

    end do

  end subroutine buildMajMap

  !!
  !! Update majorants for each region using material connections build up in buildMajMap subroutine
  !!
  subroutine updateMajorants(self, rand)
    class(transportOperatorIMC), intent(inout) :: self
    class(RNG), intent(inout)                  :: rand
    integer(shortInt)                          :: i, G, nMats
    character(100), parameter :: Here = 'updateMajorants (transportOperatorIMC_class.f90)'

    ! Check that subroutine should be called
    if (.not. allocated(self % matMajs)) return

    nMats = mm_nMat()
    G = 1   ! Can easily be extended to multiple groups later

    ! First, update array of local opacities
    do i = 1, nMats
      ! Get and verify material pointer for material i
      self % mat => mgIMCMaterial_CptrCast(self % xsData % getMaterial(i))
      if(.not.associated(self % mat)) call fatalError(Here, "Failed to get MG IMC Material")

      ! Store opacity
      self % sigmaLocal(i) = self % mat % getTotalXS(G, rand)
    end do

    ! Now update majorants for each material
    do i = 1, nMats
      self % matMajs(i) = maxval(self % sigmaLocal * self % matConnections(i, 1:nMats))
    end do

    !print *, 'Local opacities:'
    !print *, self % sigmaLocal

    !print *, 'New majorants:'
    !print *, self % matMajs

  end subroutine updateMajorants

  !!
  !! Sample position for buildMajMap subroutine (see above)
  !! Attach only necessary properties to particle:
  !!
  !!   - Position, sampled uniformly within geometry
  !!   - Direction, uniformly from unit sphere
  !!   - Type = P_PHOTON
  !!   - Group = 1, can easily extend to work with multiple groups another time
  !!
  subroutine simpleParticle(self, p, rand)
    class(transportOperatorIMC), intent(inout) :: self
    type(particle), intent(inout)              :: p
    class(RNG), intent(inout)                  :: rand
    real(defReal)                              :: mu, phi
    real(defReal), dimension(3)                :: r, dir, rand3
    integer(shortInt)                          :: matIdx, uniqueID, loops

    ! Sample points randomly within geometry until valid material is found
    loops = 0
    positionSample:do
      ! Protect against infinite loop
      loops = loops + 1
      if (loops >= 500) call fatalError(Here, '500 particles sampled in void or outside geometry')

      ! Sample position
      rand3(1) = rand % get()
      rand3(2) = rand % get()
      rand3(3) = rand % get()
      r = (self % top - self % bottom) * rand3 + self % bottom

      ! Find material under position
      call self % geom % whatIsAt(matIdx, uniqueID, r)

      ! Reject if there is no material
      if (matIdx == VOID_MAT .or. matIdx == OUTSIDE_MAT) cycle positionSample

      call p % coords % assignPosition(r)
      exit positionSample

    end do positionSample

    ! Sample Direction - chosen uniformly inside unit sphere
    mu = 2 * rand % get() - 1
    phi = rand % get() * 2*pi
    dir(1) = mu
    dir(2) = sqrt(1-mu**2) * cos(phi)
    dir(3) = sqrt(1-mu**2) * sin(phi)
    call p % coords % assignDirection(dir)

    p % type = P_PHOTON
    p % G    = 1
    p % isMG = .true.

    call self % geom % placeCoord(p % coords)

  end subroutine simpleParticle

  !!
  !! Get transport settings
  !!
  !! Sample dictionary input:
  !!
  !! transportOperator {
  !!   type    transportOperatorIMC;
  !!   cutoff  0.5;
  !!   majMap  {
  !!     nParticles 500;
  !!     steps      10;
  !!   }
  !! }
  !!
  !! As an alternative to 'steps' can specify 'lengthScale' and then steps is calculated
  !! automatically as steps = c*dt/lengthScale
  !!
  subroutine init(self, dict, geom)
    class(transportOperatorIMC), intent(inout)     :: self
    class(dictionary), intent(in)                  :: dict
    class(geometry), pointer, intent(in), optional :: geom
    class(dictionary), pointer                     :: tempDict
    integer(shortInt)                              :: nMats
    real(defReal), dimension(6)                    :: bounds
    real(defReal)                                  :: lengthScale

    ! Initialise superclass
    call init_super(self, dict)

    self % geom => geom

    ! Get timestep size
    call dict % get(self % deltaT, 'deltaT')

    ! Get cutoff value
    call dict % getOrDefault(self % cutoff, 'cutoff', 0.3_defReal)

    ! Preparation for majorant reduction subroutine
    if (dict % isPresent('majMap')) then

      if (self % cutoff == 0) call fatalError(Here, 'No need to use majorant map &
                                                    &without delta tracking')

      ! Get settings
      tempDict => dict % getDictPtr('majMap')
      call tempDict % get(self % majMapN, 'nParticles')

      if (tempDict % isPresent('steps')) then
        call tempDict % get(self % steps, 'steps')
      else
        call tempDict % get(lengthScale, 'lengthScale')
        self % steps = ceiling(lightSpeed*self % deltaT/lengthScale)
      end if

      nMats = mm_nMat()

      ! Allocate arrays
      allocate(self % matMajs(nMats))
      allocate(self % sigmaLocal(nMats))
      allocate(self % matConnections(nMats, nMats))
      self % matMajs    = 0
      self % sigmaLocal = 0

      ! Set bounding region for particle sourcing
      bounds = self % geom % bounds()
      self % bottom = bounds(1:3)
      self % top    = bounds(4:6)
    end if

  end subroutine init

end module transportOperatorIMC_class
