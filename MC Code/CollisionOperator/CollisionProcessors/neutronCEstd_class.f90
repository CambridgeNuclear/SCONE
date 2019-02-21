module neutronCEstd_class

  use numPrecision
  use endfConstants
  use genericProcedures,             only : fatalError, rotateVector, numToChar
  use dictionary_class,              only : dictionary
  use RNG_class,                     only : RNG

  ! Particle types
  use particle_class,                only : particle, phaseCoord, printType, P_NEUTRON
  use particleDungeon_class,         only : particleDungeon
  use collisionProcessor_inter,      only : collisionProcessor, collisionData ,init_super => init

  ! Nuclear Data
  use nuclearData_inter,             only : nuclearData
  use perNuclideNuclearDataCE_inter, only : perNuclideNuclearDataCE

  ! Cross-section packages to interface with nuclear data
  use xsNucMacroSet_class,    only : xsNucMacroSet_ptr
  use xsMainSet_class,        only : xsMainSet_ptr
  use xsMacroSet_class,       only : xsMacroSet_ptr
  use scatteringKernels_func, only : asymptoticScatter, targetVelocity_constXS, &
                                     asymptoticInelasticScatter
  implicit none
  private

  !!
  !! Standard (default) collision processor for CE neutrons
  !!   -> Preforms implicit fission site generation
  !!   -> Preforms analog capture
  !!   -> Treats fission as capture (only implicit generation of 2nd-ary neutrons)
  !!   -> Does not create secondary non-neutron projectiles
  !!
  !! Settings:
  !!  minE    -> minimum energy cut-off [MeV] (default = 1.0E-11)
  !!  maxE    -> maximum energy. Higher energies are set to maximum (not re-rolled) [MeV]
  !!             (default = 20.0)
  !!  tresh_E -> Energy treshold for explicit treatment of target nuclide movement [-].
  !!             Target movment is sampled if neutron energy E < kT * tresh_E where
  !!             kT is target material temperature in [MeV]. (default = 400.0)
  !!  tresh_A -> Mass treshold for explicit tratment of target nuclide movement [Mn].
  !!             Target movment is sampled if target mass A < tresh_A. (default = 1.0)
  !!
  !! Sample dictionary input:
  !!   collProcName {
  !!   #minEnergy <real>;#
  !!   #maxEnergy <real>;#
  !!   #energyTreshold <real>;#
  !!   #massTreshold <real>;~
  !!   }
  !!
  type, public, extends(collisionProcessor) :: neutronCEstd
    private
    class(perNuclideNuclearDataCE), pointer :: xsData => null() !! Nuclear Data block pointer
    !! Settings
    real(defReal) :: minE
    real(defReal) :: maxE
    real(defReal) :: tresh_E
    real(defReal) :: tresh_A

  contains
    ! Initialisation procedure
    procedure :: init

    ! Implementation of customisable procedures
    procedure :: sampleCollision
    procedure :: implicit
    procedure :: scatter
    procedure :: capture
    procedure :: fission
    procedure :: cutoffs

    ! Local procedures
    procedure,private :: elastic
    procedure,private :: inelastic
    procedure,private :: N_XN
    procedure,private :: scatterFromFixed
    procedure,private :: scatterFromMoving
    procedure,private :: scatterInLAB
  end type neutronCEstd

contains

  !!
  !! Initialise from dictionary
  !!
  subroutine init(self, dict)
    class(neutronCEstd), intent(inout) :: self
    class(dictionary), intent(in)      :: dict
    character(100), parameter :: Here = 'init (neutronCEstd_class.f90)'

    ! Call superclass
    call init_super(self, dict)

    ! Read settings for neutronCEstd
    ! Maximum and minimum energy
    call dict % getOrDefault(self % minE,'minEnergy',1.0E-11_defReal)
    call dict % getOrDefault(self % maxE,'maxEnergy',20.0_defReal)

    ! Thermal scattering kernel thresholds
    call dict % getOrDefault(self % tresh_E, 'energyTreshold', 400.0_defReal)
    call dict % getOrDefault(self % tresh_A, 'massTreshold', 1.0_defReal)

    ! Verify settings
    if( self % minE < ZERO ) call fatalError(Here,'-ve minEnergy')
    if( self % maxE < ZERO ) call fatalError(Here,'-ve maxEnergy')
    if( self % minE >= self % maxE) call fatalError(Here,'minEnergy >= maxEnergy')
    if( self % tresh_E < 0) call fatalError(Here,' -ve energyTreshold')
    if( self % tresh_A < 0) call fatalError(Here,' -ve massTreshold')

  end subroutine init

  !!
  !! Samples collision without any implicit treatment
  !!
  subroutine sampleCollision(self, p, collDat, thisCycle, nextCycle)
    class(neutronCEstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle
    type(xsNucMacroSet_ptr)              :: nucXSs
    type(xsMainSet_ptr)                  :: microXss
    real(defReal)                        :: r
    character(100),parameter :: Here = 'sampleCollision (perNuclideCollisionOpCE_class.f90)'

    ! Verify that particle is MG neutron
    if( p % isMG .or. p % type /= P_NEUTRON) then
      call fatalError(Here, 'Supports only MG Neutron. Was given CE '//printType(p % type))
    end if

    ! Verify and load nuclear data pointer
    associate ( xs => p % xsData)
      select type(xs)
        class is(perNuclideNuclearDataCE)
          self % xsData => xs

        class default
          call fatalError(Here, 'Unsupported type of Nuclear Data interface. &
                               & Only perNuclideNuclearDataCE is accepted')
      end select
    end associate

    ! Select collision nuclide
    call self % xsData % getNucMacroXs(nucXSs, p % E, collDat % matIdx)

    r = p % pRNG % get()
    collDat % nucIdx = nucXSs % invert(r)

    ! Select Main reaction channel
    call self % xsData % getMainNucXs(microXss, p % E, collDat % nucIdx)

    r = p % pRNG % get()
    collDat % MT = microXss % invert(r)

  end subroutine sampleCollision

  !!
  !! Perform implicit treatment
  !!
  subroutine implicit(self, p, collDat, thisCycle, nextCycle)
    class(neutronCEstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle
    type(xsMainSet_ptr)                  :: nuclideXss
    type(phaseCoord)                     :: pTemp
    real(defReal),dimension(3)           :: r, dir
    integer(shortInt)                    :: n, i
    real(defReal)                        :: nu, wgt, w0, rand1, E_out, mu, phi
    real(defReal)                        :: sig_fiss, sig_tot, k_eff

    ! Generate fission sites if nuclide is fissile
    if ( self % xsData % isFissileNuc(collDat % nucIdx)) then

      ! Obtain required data
      wgt   = p % w                ! Current weight
      w0    = p % preHistory % wgt ! Starting weight
      k_eff = p % k_eff            ! k_eff for normalisation
      rand1 = p % pRNG % get()     ! Random number to sample sites

      nu    = self % xsData % releaseAt(p % E, N_fission, collDat % nucIdx)
      call self % xsData % getMainNucXS(nuclideXss, p % E, collDat % nucIdx)

      sig_fiss = nuclideXss % fission()
      sig_tot  = nuclideXss % total()

      r   = p % rGlobal()
      dir = p % dirGlobal()

      ! Sample number of fission sites generated
      ! Support -ve weight particles
      n = int(abs( (wgt * nu * sig_fiss) / (w0 * sig_tot * k_eff) + rand1), shortInt)
      !n = int(abs(wgt/w0) * nu * sig_fiss/(sig_tot*k_eff) + rand1, shortInt)

      ! Store new sites in the next cycle dungeon
      wgt =  sign(w0, wgt)

      do i=1,n
          call self % xsData % sampleMuEout(mu, E_out, p % E, p % pRNG, N_fission, collDat % nucIdx)
          phi = TWO*PI * p % pRNG % get()
          dir = rotateVector(dir, mu, phi)

          if (E_out > self % maxE) E_out = self % maxE
          ! Copy extra detail from parent particle (i.e. time, flags ect.)
          pTemp       = p

          ! Overwrite position, direction, energy and weight
          pTemp % r   = r
          pTemp % dir = dir
          pTemp % E   = E_out
          pTemp % wgt = wgt

          call nextCycle % detain(pTemp)
      end do
    end if

  end subroutine implicit

  !!
  !! Process scattering reaction
  !!
  subroutine scatter(self, p, collDat, thisCycle, nextCycle)
    class(neutronCEstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle
    real(defReal)                        :: r
    character(100), parameter   :: Here ='scatter (perNuclideCollisionOp_class.f90)'

    ! Sample MT of scattering reaction. Replace lumped MT already in collDat
    r = p % pRNG % get()
    collDat % MT = self % xsData % invertScattering(p % E, r, collDat % nucIdx )

    select case(collDat % MT)
      case(N_N_elastic)
        call self % elastic(p, collDat, thisCycle, nextCycle)

      case(N_Nl1:N_Nl40, N_Ncont, N_Na, N_Np)
        call self % inelastic(p, collDat, thisCycle, nextCycle)

      case(N_2N, N_3N, N_4N)
        call self % N_XN(p, collDat, thisCycle, nextCycle)

      case default
        call fatalError(Here,'Unrecognised scattering MT number: '//numToChar(collDat % MT))

      end select

  end subroutine scatter

  !!
  !! Process capture reaction
  !!
  subroutine capture(self, p, collDat, thisCycle, nextCycle)
    class(neutronCEstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle

    p % isDead =.true.

  end subroutine capture

  !!
  !! Process fission reaction
  !!
  subroutine fission(self, p, collDat, thisCycle, nextCycle)
    class(neutronCEstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle

    p % isDead =.true.

  end subroutine fission

  !!
  !! Process elastic scattering
  !!
  subroutine elastic(self, p, collDat, thisCycle, nextCycle)
    class(neutronCEstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle
    integer(shortInt)                    :: MT, nucIdx

    ! Copy data for clarity
    MT     = collDat % MT
    nucIdx = collDat % nucIdx

    ! Scatter particle
    if (self % xsData % isInCMFrame(MT, nucIdx)) then
      collDat % A =  self % xsData % getMass(nucIdx)
      collDat % kT = self % xsData % getkT(nucIdx)

      ! Apply criterion for Free-Gas vs Fixed Target scattering
      if ((p % E > collDat % kT * self % tresh_E) .and. (collDat % A > self % tresh_A)) then
        call self % scatterFromFixed(p, collDat)

      else
        call self % scatterFromMoving(p, collDat)

      end if

    else
      call self % scatterInLAB(p, collDat)

    end if

  end subroutine elastic

  !!
  !! Process inelastic scattering
  !!
  subroutine inelastic(self, p, collDat, thisCycle, nextCycle)
    class(neutronCEstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle
    integer(shortInt)                    :: MT, nucIdx
    character(100),parameter  :: Here =' inelastic (neutronCEstd_class.f90)'

    ! Copy data for clarity
    MT     = collDat % MT
    nucIdx = collDat % nucIdx

    ! Scatter particle
    if (self % xsData % isInCMFrame(MT, nucIdx)) then
      collDat % A =  self % xsData % getMass(nucIdx)
      call self % scatterFromFixed(p, collDat)

    else
      call self % scatterInLAB(p, collDat)

    end if

  end subroutine inelastic

  !!
  !! Process N_XN scattering
  !!
  subroutine N_XN(self, p, collDat, thisCycle, nextCycle)
    class(neutronCEstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle
    integer(shortInt)                    :: MT, nucIdx
    character(100),parameter  :: Here =' N_XN (neutronCEstd_class.f90)'

    ! Copy data for clarity
    MT     = collDat % MT
    nucIdx = collDat % nucIdx

    ! Scatter particle
    if (self % xsData % isInCMFrame(MT, nucIdx)) then
      collDat % A =  self % xsData % getMass(nucIdx)
      call self % scatterFromFixed(p, collDat)

    else
      call self % scatterInLAB(p, collDat)

    end if

    ! Change particle weight
    select case(MT)
      case(N_2N)
        p % w = p % w * 2.0_defReal

      case(N_3N)
        p % w = p % w * 3.0_defReal

      case(N_4N)
        p % w = p % w * 4.0_defReal

      case default
        call fatalError(Here,'Unknown N_XN scattering. WTF?')

    end select
  end subroutine N_XN

  !!
  !! Apply cutoffs
  !!
  subroutine cutoffs(self, p, collDat, thisCycle, nextCycle)
    class(neutronCEstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle

    if (p % E < self % minE ) p % isDead = .true.

  end subroutine cutoffs

  !!
  !! Subroutine to perform scattering in LAB frame
  !! Returns mu -> cos of deflection angle in LAB frame
  !!
  subroutine scatterInLAB(self, p, collDat)
    class(neutronCEstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    real(defReal)                        :: phi    ! Azimuthal scatter angle
    real(defReal)                        :: E_out, mu
    integer(shortInt)                    :: MT, nucIdx

    ! Read data
    MT = collDat % MT
    nucIdx = collDat % nucIdx

    ! Sample scattering angles and post-collision energy
    call self % xsData % sampleMuEout(mu, E_out, p % E, p % pRNG, MT, nucIdx)
    phi = TWO * PI * p % pRNG % get()

    ! Update neutron state
    p % E = E_out
    call p % rotate(mu, phi)
    collDat % muL = mu

  end subroutine scatterInLAB

  !!
  !! Subroutine to perform scattering from stationary target.
  !! Returns mu -> cos of deflection angle in LAB frame
  !!
  subroutine scatterFromFixed(self, p, collDat)
    class(neutronCEstd), intent(inout) :: self
    class(particle), intent(inout)     :: p
    type(collisionData), intent(inout) :: collDat
    real(defReal)                      :: phi
    real(defReal)                      :: E_out
    real(defReal)                      :: E_outCM, mu
    integer(shortInt)                  :: MT, nucIdx

    ! Read data
    MT     = collDat % MT
    nucIdx = collDat % nucIdx

    ! Sample mu and outgoing energy
    call self % xsData % sampleMuEout(mu, E_outCM, p % E, p % pRNG, MT, nucIdx)

    ! Save incident energy
    E_out = p % E

    if( MT == N_N_elastic) then
      call asymptoticScatter(E_out, mu, collDat % A)

    else
      call asymptoticInelasticScatter(E_out, mu, E_outCM, collDat % A)

    end if

    ! Sample azimuthal angle
    phi = TWO * PI * p % pRNG % get()

    ! Update particle state
    call p % rotate(mu, phi)
    p % E = E_out
    collDat % muL = mu

  end subroutine scatterFromFixed

  !!
  !! Subroutine to perform scattering from moving target
  !! Supports only elastic collisions
  !!
  subroutine scatterFromMoving(self, p, collDat)
    class(neutronCEstd), intent(inout) :: self
    class(particle), intent(inout)     :: p
    type(collisionData),intent(inout)  :: collDat
    integer(shortInt)                  :: MT, nucIdx
    real(defReal)                      :: A, kT, mu
    real(defReal),dimension(3)         :: V_n           ! Neutron velocity (vector)
    real(defReal)                      :: U_n           ! Neutron speed (scalar)
    real(defReal),dimension(3)         :: dir_pre       ! Pre-collision direction
    real(defReal),dimension(3)         :: dir_post      ! Post-collicion direction
    real(defReal),dimension(3)         :: V_t, V_cm     ! Target and CM velocity
    real(defReal)                      :: phi

    ! Read data
    MT     = collDat % MT
    nucIdx = collDat % nucIdx
    A      = collDat % A
    kT     = collDat % kT

    ! Get neutron direction and velocity
    dir_pre = p % dirGlobal()
    V_n     = dir_pre * sqrt(p % E)

    ! Sample velocity of target
    V_t = targetVelocity_constXS(p % E, dir_pre, A, kT, p % pRNG)

    ! Calculate Centre-of-Mass velocity
    V_cm = (V_n + V_t *A)/(A+1)

    ! Move Neutron velocity to CM frame, store speed and calculate new normalised direction
    V_n = V_n - V_cm
    U_n = norm2(V_n)
    V_n = V_n / U_n

    ! Sample mu and phi in CM frame
    call self % xsData % sampleMu(mu, p % E, p % pRNG, MT, nucIdx)
    phi = TWO * PI * p % pRNG % get()

    ! Obtain post collision speed
    V_n = rotateVector(V_n, mu, phi) * U_n

    ! Return to LAB frame
    V_n = V_n + V_cm

    ! Calculate new neutron speed and direction
    U_n = norm2(V_n)
    dir_post = V_n / U_n

    ! Update particle state and calculate mu in LAB frame
    p % E = U_n * U_n
    call p % point(dir_post)
    collDat % muL = dot_product(dir_pre, dir_post)

  end subroutine scatterFromMoving

    
end module neutronCEstd_class
