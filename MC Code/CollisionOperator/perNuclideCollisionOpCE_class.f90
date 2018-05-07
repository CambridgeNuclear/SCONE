module perNuclideCollisionOpCE_class

  use numPrecision
  use endfConstants

  use perNuclideNuclearDataCE_inter, only : perNuclideNuclearDataCE

  use numPrecision
  use endfConstants
  use genericProcedures,     only : fatalError, rotateVector
  use RNG_class,             only : RNG
  use particle_class,        only : particle
  use particleDungeon_class, only : particleDungeon
  use collisionOperatorBase_inter,   only : collisionOperatorBase

  ! Cross-section packages to interface with nuclear data
  use xsNucMacroSet_class,    only : xsNucMacroSet_ptr
  use xsMainSet_class,        only : xsMainSet_ptr

  use scatteringKernels_func, only : asymptoticScatter, targetVelocity_constXS

  implicit none
  private

  real(defReal), parameter :: energyCutoff = 2.0E-11, energyMaximum = 20.0


  type, public,extends(collisionOperatorBase)  :: perNuclideCollisionOpCE
   !*DEBUG private
    real(defReal)                           :: E      !! Particle energy
    class(perNuclideNuclearDataCE), pointer :: xsData !! Nuclear Data block pointer
  contains
    procedure :: sampleCollision
    procedure :: implicit
    procedure :: scatter
    procedure :: elastic
    procedure :: inelastic
    procedure :: N_XN
    procedure :: capture
    procedure :: fission
    procedure :: cutoffs

    ! * Private Procedures
    procedure,private :: scatterFromFixed
    procedure,private :: scatterFromMoving
    procedure,private :: scatterInLAB

  end type perNuclideCollisionOpCE

contains

  subroutine sampleCollision(self,p,thisCycle,nextCycle)
    class(perNuclideCollisionOpCE), intent(inout) :: self
    class(particle), intent(inout)                :: p
    class(particleDungeon),intent(inout)          :: thisCycle
    class(particleDungeon),intent(inout)          :: nextCycle
    type(xsNucMacroSet_ptr)                       :: nucXSs
    type(xsMainSet_ptr)                           :: microXss
    real(defReal)                                 :: r

    ! Load collision energy
    self % E = p % E

    ! Select collision nuclide
    call self % xsData % getNucMacroXs(nucXSs, self % E, self % matIdx)

    r = self % pRNG % get()
    self % nucIdx = nucXSs % invert(r)

    ! Select Main reaction channel
    call self % xsData % getMainNucXs(microXss, self % E, self % nucIdx)

    r = self % pRNG % get()
    self % MT = microXss % invert(r)

  end subroutine sampleCollision
    
  subroutine implicit(self,p,thisCycle,nextCycle)
    class(perNuclideCollisionOpCE), intent(inout) :: self
    class(particle), intent(inout)                :: p
    class(particleDungeon),intent(inout)          :: thisCycle
    class(particleDungeon),intent(inout)          :: nextCycle
    type(xsMainSet_ptr)                      :: nuclideXss
    type(particle)                           :: pTemp
    real(defReal),dimension(3)               :: r, dir
    integer(shortInt)                        :: n, i
    real(defReal)                            :: nu, wgt, r1, E_out, mu, phi
    real(defReal)                            :: sig_fiss, sig_tot, k_eff

    ! Generate fission sites if nuclide is fissile
    if ( self % xsData % isFissileNuc(self % nucIdx)) then


      ! Obtain required data
      wgt   = p % w
      nu    = self % xsData % releaseAt(self % E, N_fission, self % nucIdx)
      k_eff = nextCycle % k_eff
      r1    = self % pRNG % get()
      call self % xsData % getMainNucXS(nuclideXss,self % E, self % nucIdx)

      sig_fiss = nuclideXss % fission()
      sig_tot  = nuclideXss % total()

      r   = p % rGlobal()
      dir = p % dirGlobal()

      ! Sample number of fission sites generated
      n = int(wgt * nu * sig_fiss/(sig_tot*k_eff) + r1, shortInt)

      ! Throw new sites to the next cycle dungeon
      wgt = 1.0
      do i=1,n
          call self % xsData % sampleMuEout(mu, E_out, self % E, self % pRNG, N_fission, self % nucIdx)
          phi = 2*PI * self % pRNG % get()
          dir = rotateVector(dir,mu,phi)

          if (E_out > energyMaximum) E_out = energyMaximum

          call pTemp % build(r,dir,E_out,wgt)
          call nextCycle % throw(pTemp)

      end do
    end if


  end subroutine implicit

  subroutine scatter(self,p,thisCycle,nextCycle)
    class(perNuclideCollisionOpCE), intent(inout) :: self
    class(particle), intent(inout)                :: p
    class(particleDungeon),intent(inout)          :: thisCycle
    class(particleDungeon),intent(inout)          :: nextCycle
    character(100), parameter     :: Here ='scatter (perNuclideCollisionOp_class.f90)'

    !* DO not sample number yet to preserve the sequence
    self % MT = self % xsData % invertScattering(self %E, 0.5_defReal)

    select case(self % MT)
      case(N_N_elastic)
        call self % elastic(p,thisCycle,nextCycle)

      case(N_Nl1:N_Nl40, N_Ncont)
        call self % inelastic(p,thisCycle,nextCycle)

      case(N_2N,N_3N,N_4N)
        call self % N_XN(p,thisCycle,nextCycle)

      case default
        call fatalError(Here,'Unrecognised scattering MT number')

      end select

  end subroutine scatter

  subroutine elastic(self,p,thisCycle,nextCycle)
    class(perNuclideCollisionOpCE), intent(inout) :: self
    class(particle), intent(inout)                :: p
    class(particleDungeon),intent(inout)          :: thisCycle
    class(particleDungeon),intent(inout)          :: nextCycle
    integer(shortInt)                             :: MT, nucIdx
    real(defReal)                                 :: E          ! Pre-collision energy
    real(defReal)                                 :: kT, A      ! Target temperature[MeV] and mass [Mn]
    real(defReal)                                 :: muL        ! Cosine of scattering in LAB frame

    ! Copy operator variables to local copies for clarity
    MT     = self % MT
    E      = self % E
    nucIdx = self % nucIdx

    if (self % xsData % isInCMFrame(MT, nucIdx)) then
      A =  self % xsData % getMass(nucIdx)
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

    self % muL = muL

  end subroutine elastic

  subroutine inelastic(self,p,thisCycle,nextCycle)
    class(perNuclideCollisionOpCE), intent(inout) :: self
    class(particle), intent(inout)                :: p
    class(particleDungeon),intent(inout)          :: thisCycle
    class(particleDungeon),intent(inout)          :: nextCycle
    character(100),parameter              :: Here =' inelastic (perNuclideCollisionOp_class.f90)'

    call fatalError(Here,'Inelastic Scattering is not implemented')

  end subroutine inelastic

  subroutine N_XN(self,p,thisCycle,nextCycle)
    class(perNuclideCollisionOpCE), intent(inout) :: self
    class(particle), intent(inout)                :: p
    class(particleDungeon),intent(inout)          :: thisCycle
    class(particleDungeon),intent(inout)          :: nextCycle
    character(100),parameter              :: Here =' N_XN (perNuclideCollisionOp_class.f90)'

    call fatalError(Here,'N_XN Scattering is not implemented')

  end subroutine N_XN

  subroutine capture(self,p,thisCycle,nextCycle)
    class(perNuclideCollisionOpCE), intent(inout) :: self
    class(particle), intent(inout)                :: p
    class(particleDungeon),intent(inout)          :: thisCycle
    class(particleDungeon),intent(inout)          :: nextCycle

    p % isDead =.true.

  end subroutine capture

  subroutine fission(self,p,thisCycle,nextCycle)
    class(perNuclideCollisionOpCE), intent(inout) :: self
    class(particle), intent(inout)                :: p
    class(particleDungeon),intent(inout)          :: thisCycle
    class(particleDungeon),intent(inout)          :: nextCycle

    p % isDead =.true.

  end subroutine fission

  subroutine cutoffs(self,p,thisCycle,nextCycle)
    class(perNuclideCollisionOpCE), intent(inout) :: self
    class(particle), intent(inout)                :: p
    class(particleDungeon),intent(inout)          :: thisCycle
    class(particleDungeon),intent(inout)          :: nextCycle

    if (p % E < energyCutoff ) p % isDead = .true.

  end subroutine cutoffs

  !!
  !! Subroutine to perform scattering in LAB frame
  !! Returns mu -> cos of deflection angle in LAB frame
  !!
  subroutine scatterInLAB(self,mu,p,E,MT,nucIdx)
    class(perNuclideCollisionOpCE), intent(inout) :: self
    real(defReal), intent(out)              :: mu    ! Returned deflection angle cos in LAB
    class(particle), intent(inout)           :: p
    real(defReal), intent(in)               :: E      ! Neutron energy
    integer(shortInt),intent(in)            :: MT     ! Reaction MT number
    integer(shortInt),intent(in)            :: nucIdx ! Target nuclide index
    real(defReal)                           :: phi    ! Azimuthal scatter angle
    real(defReal)                           :: E_out

    ! Sample scattering angles and post-collision energy
    call self % xsData % sampleMuEout(mu, E_out, E, self % pRNG, MT, nucIdx)
    phi = 2*PI* self % pRNG % get()

    ! Update neutron state
    p % E = E_out
    call p % rotate(mu,phi)

  end subroutine scatterInLAB

  !!
  !! Subroutine to perform scattering from stationary target.
  !! Returns mu -> cos of deflection angle in LAB frame
  !!
  subroutine scatterFromFixed(self,mu,p,E,A,MT,nucIdx)
    class(perNuclideCollisionOpCE), intent(inout) :: self
    real(defReal), intent(out)              :: mu          ! Returned deflection angle cos in LAB
    class(particle), intent(inout)           :: p
    real(defReal), intent(in)               :: E           ! Neutron energy
    real(defReal), intent(in)               :: A           ! Target weight
    integer(shortInt),intent(in)            :: MT          ! Reaction MT number
    integer(shortInt),intent(in)            :: nucIdx      ! Target nuclide index
    real(defReal)                           :: phi
    real(defReal)                           :: E_out
    character(100),parameter       :: Here = 'scatterFromStationary (collisionOperator_class.f90)'

    ! Sample mu and outgoing energy
    call self % xsData % sampleMuEout(mu, E_out, E, self % pRNG, MT, nucIdx)

    select case(MT)
      case(N_N_elastic)
        call asymptoticScatter(E_out,mu,A)

      case default
        call fatalError(Here,'Unknown MT number')

    end select

    ! Sample azimuthal angle
    phi = 2*PI * self % pRNG % get()

    ! Update particle state
    call p % rotate(mu,phi)
    p % E = E_out

  end subroutine scatterFromFixed


  !!
  !! Subroutine to perform scattering  from moving target
  !! Returns mu -> cos of deflection angle in LAB
  !!
  subroutine scatterFromMoving(self,mu,p,E,A,kT,MT,nucIdx)
    class(perNuclideCollisionOpCE), intent(inout) :: self
    real(defReal), intent(out)                    :: mu
    class(particle), intent(inout)          :: p
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
    dir_pre = p % dirGlobal()
    V_n     = dir_pre * sqrt(E)

    ! Sample velocity of target
    V_t = targetVelocity_constXS(E, dir_pre, A, kT, self % pRNG)

    ! Calculate Centre-of-Mass velocity
    V_cm = (V_n + V_t *A)/(A+1)

    ! Move Neutron velocity to CM frame, store speed and calculate new normalised direction
    V_n = V_n - V_cm
    U_n = norm2(V_n)
    V_n = V_n / U_n

    ! Sample mu and phi in CM frame
    call self % xsData % sampleMu(mu, E, self % pRNG, MT, nucIdx)
    phi = 2*PI*self % pRNG % get()

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




end module perNuclideCollisionOpCE_class
