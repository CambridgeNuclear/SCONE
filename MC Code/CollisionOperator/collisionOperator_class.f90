module collisionOperator_class

  use numPrecision
  use endfConstants
  use genericProcedures,      only : fatalError, rotateVector
  use RNG_class,              only : RNG
  use particle_class,         only : particle
  use particleDungeon_class,  only : particleDungeon
  use byNucNoMT_class,        only : byNucNoMT

  ! Cross-section packages to interface with nuclear data
  use matNucCDF_class,        only : matNucCDF
  use xsMainCDF_class,        only : xsMainCDF
  use xsMainSet_class,        only : xsMainSet

  use scatteringKernels_func, only : asymptoticScatter, targetVelocity_constXS

  implicit none
  private

  !!
  !! ***** Eventual improvement
  !! *** Define a common-block like type to contain all important data for collision procesing
  !! ** Move all physic subroutines to "a library" in external module. Use the type to interface
  !! * with it. Thus collision operators will become more higher level. Tell what you want to do
  !! Do not tell how to do it!
  !!
  real(defReal), parameter :: energyCutoff = 2.0E-11, energyMaximum = 20.0

  type, public :: collisionOperator
   !* private ** DEBUG
    class(byNucNoMT), pointer :: xsData => null()
    class(RNG), pointer       :: locRNG => null()
  contains
    procedure :: attachXsData
    procedure :: collide

    ! Private procedures
    procedure :: performScattering
    procedure :: performCapture
    procedure :: performFission
    procedure :: scatterFromFixed
    procedure :: scatterFromMoving
    procedure :: scatterInLAB
    procedure :: createFissionSites

  end type collisionOperator

contains

  !!
  !! Initialise XS operator by providing pointer to XS data block
  !!
  subroutine attachXsData(self,xsData)
    class(collisionOperator), intent(inout) :: self
    class(byNucNoMT),pointer, intent(in)    :: xsData
    character(100), parameter               :: Here =' attachXsData (collisionOperator_class.f90)'

    if(.not.associated(xsData)) call fatalError(Here,'Allocated xs data must be provided')

    self % xsData => xsData

  end subroutine attachXSData


  !!
  !! Subroutine to collide a neutron. Chooses collision nuclide and main reaction channel.
  !! Calls approperiate procedure to change neutron state
  !!
  subroutine collide(self,p,thisCycle,nextCycle)
    class(collisionOperator), intent(inout) :: self
    type(particle), intent(inout)           :: p
    type(particleDungeon),intent(inout)     :: thisCycle
    type(particleDungeon),intent(inout)     :: nextCycle
    real(defReal)                           :: E
    integer(shortInt)                       :: matIdx
    integer(shortInt)                       :: nucIdx
    integer(shortInt)                       :: MT
    type(matNucCDF),pointer                 :: nuclideCDF
    type(xsMainCDF),pointer                 :: reactionCDF
    real(defReal)                           :: r

    ! Retrive Pointer to the random number generator
    self % locRNG => p % pRNG

    ! Load neutron energy and material
    E = p % E
    matIdx = p % matIdx

    ! Select collision nuclide
    call self % xsData % getMatNucCDF(nuclideCDF, E, matIdx)

    r = self % locRNG % get()
    nucIdx = nuclideCDF % invert(r)

    ! Select Main reaction channel
    call self % xsData % getMainNucCdf(reactionCDF, E, nucIdx)

    r = self % locRNG % get()
    MT = reactionCDF % invert(r)

    ! Generate fission sites if nuclide is fissile

    if ( self % xsData % isFissile(nucIdx)) then
      call self % createFissionSites(nextCycle,p,nucIdx)

    end if

    ! Call procedure to do reaction processing
    select case (MT)
      case(anyScatter)
        call self % performScattering(p,nucIdx)

      case(anyCapture)
        call self % performCapture(p,nucIdx)

      case(anyFission)
        call self % performFission(p,nucIdx)

    end select

    if (p % E < energyCutoff ) p % isDead = .true.

  end subroutine collide

  !!
  !! Change particle state in scattering reaction
  !! Critarion for Free-Gas vs Fixed Target scattering is taken directly from MCNP manual chapter 2
  !!
  subroutine performScattering(self,p,nucIdx)
    class(collisionOperator), intent(inout) :: self
    type(particle), intent(inout)           :: p
    integer(shortInt),intent(in)            :: nucIdx
    integer(shortInt)                       :: MT
    real(defReal)                           :: E          ! Pre-collision energy
    real(defReal)                           :: kT, A      ! Target temperature[MeV] and mass [Mn]
    real(defReal)                           :: muL        ! Cosine of scattering in LAB frame

    ! Assign MT number -> only elastic scattering at this moment
    MT = N_N_elastic

    E = p % E


    if (self % xsData % isInCMFrame(MT, nucIdx)) then
      A =  self % xsData % getWeight(nucIdx)
      kT = self % xsData % getkT(nucIdx)

      ! Apply criterion for Free-Gas vs Fixed Target scattering
      if ((E > kT*400.0) .and. (A>1.0)) then
        call self % scatterFromFixed(muL,p,E,A,MT,nucIdx)

      else
        call self % scatterFromMoving(muL,p,E,A,kT,MT,nucIdx)

      end if

    else
      call self % scatterInLAB(muL,p,E,MT,nucIdx)

    end if

  end subroutine performScattering


  !!
  !! Subroutine to perform scattering in LAB frame
  !! Returns mu -> cos of deflection angle in LAB frame
  !!
  subroutine scatterInLAB(self,mu,p,E,MT,nucIdx)
    class(collisionOperator), intent(inout) :: self
    real(defReal), intent(out)              :: mu    ! Returned deflection angle cos in LAB
    type(particle), intent(inout)           :: p
    real(defReal), intent(in)               :: E      ! Neutron energy
    integer(shortInt),intent(in)            :: MT     ! Reaction MT number
    integer(shortInt),intent(in)            :: nucIdx ! Target nuclide index
    real(defReal)                           :: phi    ! Azimuthal scatter angle
    real(defReal)                           :: E_out

    ! Sample scattering angles and post-collision energy
    call self % xsData % sampleMuEout(mu, E_out, E, self % locRNG, MT, nucIdx)
    phi = 2*PI* self % locRNG % get()

    ! Update neutron state
    p % E = E_out
    call p % rotate(mu,phi)

  end subroutine scatterInLAB

  !!
  !! Subroutine to perform scattering from stationary target.
  !! Returns mu -> cos of deflection angle in LAB frame
  !!
  subroutine scatterFromFixed(self,mu,p,E,A,MT,nucIdx)
    class(collisionOperator), intent(inout) :: self
    real(defReal), intent(out)              :: mu          ! Returned deflection angle cos in LAB
    type(particle), intent(inout)           :: p
    real(defReal), intent(in)               :: E           ! Neutron energy
    real(defReal), intent(in)               :: A           ! Target weight
    integer(shortInt),intent(in)            :: MT          ! Reaction MT number
    integer(shortInt),intent(in)            :: nucIdx      ! Target nuclide index
    real(defReal)                           :: phi
    real(defReal)                           :: E_out
    character(100),parameter       :: Here = 'scatterFromStationary (collisionOperator_class.f90)'

    ! Sample mu and outgoing energy
    call self % xsData % sampleMuEout(mu, E_out, E, self % locRNG, MT, nucIdx)

    select case(MT)
      case(N_N_elastic)
        call asymptoticScatter(E_out,mu,A)

      case default
        call fatalError(Here,'Unknown MT number')

    end select

    ! Sample azimuthal angle
    phi = 2*PI * self % locRNG % get()

    ! Update particle state
    call p % rotate(mu,phi)
    p % E = E_out

  end subroutine scatterFromFixed


  !!
  !! Subroutine to perform scattering  from moving target
  !! Returns mu -> cos of deflection angle in LAB
  !!
  subroutine scatterFromMoving(self,mu,p,E,A,kT,MT,nucIdx)
    class(collisionOperator), intent(inout) :: self
    real(defReal), intent(out)              :: mu
    type(particle), intent(inout)           :: p
    real(defReal), intent(in)               :: E            ! Neutron energy
    real(defReal), intent(in)               :: A
    real(defReal),intent(in)                :: kT           ! Target temperature
    integer(shortInt),intent(in)            :: MT
    integer(shortInt),intent(in)            :: nucIdx        ! Target nuclide index
    real(defReal),dimension(3)              :: V_n           ! Neutron velocity (vector)
    real(defReal)                           :: U_n           ! Neutron speed (scalar)
    real(defReal),dimension(3)              :: dir_pre       ! Pre-collision direction
    real(defReal),dimension(3)              :: dir_post      ! Post-collicion direction
    real(defReal),dimension(3)              :: V_t, V_cm     ! Target and CM velocity
    real(defReal)                           :: phi

    ! Get neutron direction and velocity
    dir_pre = p % globalDir()
    V_n     = dir_pre * sqrt(E)

    ! Sample velocity of target
    V_t = targetVelocity_constXS(E, dir_pre, A, kT, self % locRNG)

    ! Calculate Centre-of-Mass velocity
    V_cm = (V_n + V_t *A)/(A+1)

    ! Move Neutron velocity to CM frame, store speed and calculate new normalised direction
    V_n = V_n - V_cm
    U_n = norm2(V_n)
    V_n = V_n / U_n

    ! Sample mu and phi in CM frame
    call self % xsData % sampleMu(mu, E, self % locRNG, MT, nucIdx)
    phi = 2*PI*self % locRNG % get()

    ! Obtain post collision speed
    V_n = rotateVector(V_n,mu,phi) * U_n

    ! Return to LAB frame
    V_n = V_n + V_cm

    ! Calculate new neutron speed and direction
    U_n = norm2(V_n)
    dir_post = V_n / U_n

    ! Update particle state and calculate mu in LAB frame
    p % E = U_n * U_n
    call p % point(dir_post)
    mu = dot_product(dir_pre,dir_post)

  end subroutine scatterFromMoving




  !!
  !! Change particle state in capture reaction
  !!
  subroutine performCapture(self,p,nucIdx)
    class(collisionOperator), intent(inout) :: self
    type(particle), intent(inout)           :: p
    integer(shortInt),intent(in)            :: nucIdx

    p % isDead =.true.

  end subroutine performCapture

  !!
  !! Change particle state in fission reaction
  !!
  subroutine performFission(self,p,nucIdx)
    class(collisionOperator), intent(inout) :: self
    type(particle), intent(inout)           :: p
    integer(shortInt),intent(in)            :: nucIdx

    p % isDead =.true.

  end subroutine performFission

  !!
  !!
  !!
  subroutine createFissionSites(self,nextCycle,p,nucIdx)
    class(collisionOperator), intent(inout)  :: self
    type(particleDungeon), intent(inout)     :: nextCycle
    type(particle), intent(in)               :: p
    integer(shortInt), intent(in)            :: nucIdx
    type(xsMainSet),pointer                  :: nuclideXss
    type(particle)                           :: pTemp
    real(defReal),dimension(3)               :: r, dir
    integer(shortInt)                        :: n, i
    real(defReal)                            :: E, nu, wgt, r1, E_out, mu, phi
    real(defReal)                            :: sig_fiss, sig_tot, k_eff

    ! Obtain required data
    E     = p % E
    wgt   = p % w
    nu    = self % xsData % releaseAt(E, N_fission, nucIdx)
    k_eff = nextCycle % k_eff
    r1    = self % locRNG % get()
    call self % xsData % getMainNucXS(nuclideXss,E,nucIdx)

    sig_fiss = nuclideXss % fission
    sig_tot  = nuclideXss % total

    r   = p % rGlobal()
    dir = p % globalDir()

    ! Sample number of fission sites generated
    n = int(wgt * nu * sig_fiss/(sig_tot*k_eff) + r1, shortInt)

    ! Throw new sites to the next cycle dungeon
    wgt = 1.0
    do i=1,n
      call self % xsData % sampleMuEout(mu, E_out, E, self % locRNG, N_fission, nucIdx)
      phi = 2*PI * self % locRNG % get()
      dir = rotateVector(dir,mu,phi)

      if (E_out > energyMaximum) E_out = energyMaximum

      call pTemp % build(r,dir,E_out,wgt)
      call nextCycle % throw(pTemp)
    end do


  end subroutine createFissionSites



end module collisionOperator_class
