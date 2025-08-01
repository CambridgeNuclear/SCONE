module neutronCEstd_class

  use numPrecision
  use endfConstants
  use universalVariables,            only : REJECTED, INF
  use genericProcedures,             only : fatalError, rotateVector, numToChar
  use dictionary_class,              only : dictionary
  use RNG_class,                     only : RNG

  ! Particle types
  use particle_class,                only : particle, particleState, printType, P_NEUTRON, P_PRECURSOR
  use particleDungeon_class,         only : particleDungeon

  ! Abstarct interface
  use collisionProcessor_inter,      only : collisionProcessor, collisionData ,init_super => init

  ! Nuclear Data Interfaces
  use nuclearDataReg_mod,            only : ndReg_getNeutronCE => getNeutronCE
  use nuclearDatabase_inter,         only : nuclearDatabase
  use ceNeutronDatabase_inter,       only : ceNeutronDatabase
  use ceNeutronMaterial_class,       only : ceNeutronMaterial, ceNeutronMaterial_CptrCast
  use ceNeutronNuclide_inter,        only : ceNeutronNuclide, ceNeutronNuclide_CptrCast

  ! Nuclear reactions
  use reactionHandle_inter,          only : reactionHandle
  use uncorrelatedReactionCE_inter,  only : uncorrelatedReactionCE, uncorrelatedReactionCE_CptrCast
  use neutronScatter_class,          only : neutronScatter, neutronScatter_TptrCast
  use fissionCE_class,               only : fissionCE, fissionCE_TptrCast

  ! Cross-Section Packages
  use neutronXsPackages_class,       only : neutronMicroXSs

  ! Scattering procedures
  use scatteringKernels_func, only : asymptoticScatter, targetVelocity_constXS, &
                                     asymptoticInelasticScatter, targetVelocity_DBRCXS, &
                                     relativeEnergy_constXS

  ! Tally interfaces
  use tallyAdmin_class,       only : tallyAdmin

  implicit none
  private

  !!
  !! Standard (default) scalar collision processor for CE neutrons
  !!   -> Preforms implicit fission site generation
  !!   -> Preforms analog capture
  !!   -> Treats fission as capture (only implicit generation of 2nd-ary neutrons)
  !!   -> Does not create secondary non-neutron projectiles
  !!
  !! Settings:
  !!  minE       -> minimum energy cut-off [MeV] (default = 1.0E-11)
  !!  maxE       -> maximum energy. Higher energies are set to maximum (not re-rolled) [MeV]
  !!                (default = 20.0)
  !!  threshE    -> Energy threshold for explicit treatment of target nuclide movement [-].
  !!                Target movement is sampled if neutron energy E < kT * threshE where
  !!                kT is target material temperature in [MeV]. (default = 400.0)
  !!  threshA    -> Mass threshold for explicit treatment of target nuclide movement [Mn].
  !!               Target movement is sampled if target mass A < threshA. (default = 1.0)
  !!  DBRCeMin   -> Minimum energy to which DBRC is applied
  !!  DBRCeMax   -> Maximum energy to which DBRC is applied
  !!  makePrec   -> Produce precursor particles, used in dynamic calculations (default = false)
  !!  promptOnly -> If true, prevents delayed neutrons or precursors from being produced
  !!
  !! Sample dictionary input:
  !!   collProcName {
  !!   type             neutronCEstd;
  !!   #minEnergy       <real>;#
  !!   #maxEnergy       <real>;#
  !!   #energyThreshold <real>;#
  !!   #massThreshold   <real>;#
  !!   #makePrec        <bool>;#
  !!   #promptOnly      <bool>;#
  !!   }
  !!
  type, public, extends(collisionProcessor) :: neutronCEstd
    private
    !! Nuclear Data block pointer -> public so it can be used by subclasses (protected member)
    class(ceNeutronDatabase), pointer, public :: xsData => null()
    class(ceNeutronMaterial), pointer, public :: mat    => null()
    class(ceNeutronNuclide),  pointer, public :: nuc    => null()

    !! Settings - private
    real(defReal) :: minE
    real(defReal) :: maxE
    real(defReal) :: threshE
    real(defReal) :: threshA
    real(defReal) :: DBRCeMin
    real(defReal) :: DBRCeMax
    logical(defBool) :: makePrec = .false.
    logical(defBool) :: promptOnly = .false.

  contains
    ! Initialisation procedure
    procedure :: init

    ! Implementation of customisable procedures
    procedure :: sampleCollision
    procedure :: implicit
    procedure :: elastic
    procedure :: inelastic
    procedure :: capture
    procedure :: fission
    procedure :: cutoffs

    ! Local procedures
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
    call dict % getOrDefault(self % threshE, 'energyThreshold', 400.0_defReal)
    call dict % getOrDefault(self % threshA, 'massThreshold', 1.0_defReal)
    
    ! Precursor settings
    call dict % getOrDefault(self % makePrec, 'makePrec', .false.)
    call dict % getOrDefault(self % promptOnly, 'promptOnly', .false.)

    ! Verify settings
    if (self % minE < ZERO) call fatalError(Here,'-ve minEnergy')
    if (self % maxE < ZERO) call fatalError(Here,'-ve maxEnergy')
    if (self % minE >= self % maxE) call fatalError(Here,'minEnergy >= maxEnergy')
    if (self % threshE < 0) call fatalError(Here,' -ve energyThreshold')
    if (self % threshA < 0) call fatalError(Here,' -ve massThreshold')
    if (self % makePrec .and. self % promptOnly) call fatalError(Here,&
            'Incompatible options: cannot makePrecursors and have promptOnly neutrons!')

    ! DBRC energy limits
    call dict % getOrDefault(self % DBRCeMin,'DBRCeMin', (1.0E-8_defReal))
    call dict % getOrDefault(self % DBRCeMax,'DBRCeMax', (200E-6_defReal))

  end subroutine init

  !!
  !! Samples collision without any implicit treatment
  !!
  subroutine sampleCollision(self, p, tally, collDat, thisCycle, nextCycle)
    class(neutronCEstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(tallyAdmin), intent(inout)      :: tally
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle
    type(neutronMicroXSs)                :: microXSs
    real(defReal)                        :: r
    character(100),parameter :: Here = 'sampleCollision (neutronCEstd_class.f90)'

    ! Verify that particle is CE neutron
    if (p % isMG .or. p % type /= P_NEUTRON) then
      call fatalError(Here, 'Supports only CE Neutron. Was given MG '//printType(p % type))
    end if

    ! Verify and load nuclear data pointer
    self % xsData => ndReg_getNeutronCE()
    if (.not.associated(self % xsData)) call fatalError(Here, 'There is no active Neutron CE data!')

    ! Verify and load material pointer
    self % mat => ceNeutronMaterial_CptrCast(self % xsData % getMaterial(p % matIdx()))
    if (.not.associated(self % mat)) call fatalError(Here, 'Material is not ceNeutronMaterial')

    ! Select collision nuclide
    call self % mat % sampleNuclide(p % E, p % pRNG, collDat % nucIdx, collDat % E)

    ! If nuclide was rejected in TMS loop return to tracking
    if (collDat % nucIdx == REJECTED) then
      collDat % MT = noInteraction
      return
    end if

    self % nuc => ceNeutronNuclide_CptrCast(self % xsData % getNuclide(collDat % nucIdx))
    if (.not.associated(self % nuc)) call fatalError(Here, 'Failed to retrieve CE Neutron Nuclide')

    ! Select Main reaction channel
    call self % nuc % getMicroXSs(microXss, collDat % E, self % mat % kT, p % pRNG)
    r = p % pRNG % get()
    collDat % MT = microXss % invert(r)

  end subroutine sampleCollision

  !!
  !! Perform implicit treatment
  !!
  subroutine implicit(self, p, tally, collDat, thisCycle, nextCycle)
    class(neutronCEstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(tallyAdmin), intent(inout)      :: tally
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle
    type(fissionCE), pointer             :: fission
    type(neutronMicroXSs)                :: microXSs
    type(particleState)                  :: pTemp
    real(defReal),dimension(3)           :: r, dir
    integer(shortInt)                    :: n, i
    real(defReal)                        :: wgt, w0, rand1, E_out, mu, phi
    real(defReal)                        :: sig_nufiss, sig_tot, k_eff, lambda
    character(100),parameter             :: Here = 'implicit (neutronCEstd_class.f90)'

    ! Generate fission sites if nuclide is fissile
    if (self % nuc % isFissile()) then

      ! Obtain required data
      wgt   = p % w                ! Current weight
      w0    = p % preHistory % wgt ! Starting weight
      k_eff = p % k_eff            ! k_eff for normalisation
      rand1 = p % pRNG % get()     ! Random number to sample sites

      ! Retrieve cross section at the energy used for reaction sampling
      call self % nuc % getMicroXSs(microXSs, collDat % E, self % mat % kT, p % pRNG)

      sig_nufiss = microXSs % nuFission
      sig_tot    = microXSs % total

      ! Sample number of fission sites generated
      ! Support -ve weight particles
      n = int(abs( (wgt * sig_nufiss) / (w0 * sig_tot * k_eff)) + rand1, shortInt)

      ! Shortcut particle generation if no particles were sampled
      if (n < 1) return

      ! Get fission reaction
      fission => fissionCE_TptrCast(self % xsData % getReaction(N_FISSION, collDat % nucIdx))
      if(.not.associated(fission)) call fatalError(Here, "Failed to get fissionCE")

      ! Store new sites in the next cycle dungeon
      wgt = sign(w0, wgt)
      r   = p % rGlobal()

      do i = 1, n
        call fission % sampleOut(mu, phi, E_out, p % E, p % pRNG, lambda)
        
        ! Skip if a delayed particle is produced in prompt-only mode
        if (self % promptOnly .and. lambda < huge(lambda)) cycle
        
        dir = rotateVector(p % dirGlobal(), mu, phi)
        
        if (E_out > self % maxE) E_out = self % maxE

        ! Copy extra detail from parent particle (i.e. time, flags ect.)
        pTemp       = p

        ! Overwrite position, direction, energy and weight
        pTemp % r   = r
        pTemp % dir = dir
        pTemp % E   = E_out
        pTemp % wgt = wgt
        pTemp % collisionN = 0

        ! If storing precursors, do so when a finite lambda occurs
        if (self % makePrec .and. lambda < huge(lambda)) then
          pTemp % lambda = lambda
          pTemp % type = P_PRECURSOR

        end if

        call nextCycle % detain(pTemp)

        ! Report birth of new particle
        call tally % reportSpawn(N_FISSION, p, pTemp)

      end do

    end if

  end subroutine implicit

  !!
  !! Process capture reaction
  !!
  subroutine capture(self, p, tally, collDat, thisCycle, nextCycle)
    class(neutronCEstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(tallyAdmin), intent(inout)      :: tally
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle

    p % isDead =.true.

  end subroutine capture

  !!
  !! Process fission reaction
  !!
  subroutine fission(self, p, tally, collDat, thisCycle, nextCycle)
    class(neutronCEstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(tallyAdmin), intent(inout)      :: tally
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle

    p % isDead =.true.

  end subroutine fission

  !!
  !! Process elastic scattering
  !!
  !! All CE elastic scattering happens in the CM frame
  !!
  subroutine elastic(self, p, tally, collDat, thisCycle, nextCycle)
    class(neutronCEstd), intent(inout)     :: self
    class(particle), intent(inout)         :: p
    type(tallyAdmin), intent(inout)        :: tally
    type(collisionData), intent(inout)     :: collDat
    class(particleDungeon),intent(inout)   :: thisCycle
    class(particleDungeon),intent(inout)   :: nextCycle
    class(uncorrelatedReactionCE), pointer :: reac
    logical(defBool)                       :: isFixed, hasDBRC
    character(100),parameter :: Here = 'elastic (neutronCEstd_class.f90)'

    ! Assess if thermal scattering data is needed or not
    if (self % nuc % needsSabEl(p % E)) collDat % MT = N_N_ThermEL

    ! Get reaction
    reac => uncorrelatedReactionCE_CptrCast( self % xsData % getReaction(collDat % MT, collDat % nucIdx))
    if (.not.associated(reac)) call fatalError(Here,'Failed to get elastic neutron scatter')

    ! Scatter particle
    collDat % A =  self % nuc % getMass()

    ! Retrieve kT from either material or nuclide
    if (self % mat % useTMS(p % E)) then
      collDat % kT = self % mat % kT
    else
      collDat % kT = self % nuc % getkT()
    end if

    ! Check is DBRC is on
    hasDBRC = self % nuc % hasDBRC()

    isFixed = (.not. hasDBRC) .and. (p % E > collDat % kT * self % threshE) &
              & .and. (collDat % A > self % threshA)

    ! Apply criterion for Free-Gas vs Fixed Target scattering
    if (.not. reac % inCMFrame()) then
      call self % scatterInLAB(p, collDat, reac)
    elseif (isFixed) then
      call self % scatterFromFixed(p, collDat, reac)
    else
      call self % scatterFromMoving(p, collDat, reac)
    end if

  end subroutine elastic

  !!
  !! Process inelastic scattering
  !!
  subroutine inelastic(self, p, tally, collDat, thisCycle, nextCycle)
    class(neutronCEstd), intent(inout)     :: self
    class(particle), intent(inout)         :: p
    type(tallyAdmin), intent(inout)        :: tally
    type(collisionData), intent(inout)     :: collDat
    class(particleDungeon),intent(inout)   :: thisCycle
    class(particleDungeon),intent(inout)   :: nextCycle
    class(uncorrelatedReactionCE), pointer :: reac
    character(100),parameter  :: Here =' inelastic (neutronCEstd_class.f90)'

    ! Invert inelastic scattering and get reaction
    collDat % MT = self % nuc % invertInelastic(collDat % E, p % pRNG)
    reac => uncorrelatedReactionCE_CptrCast(self % xsData % getReaction(collDat % MT, collDat % nucIdx))
    if (.not.associated(reac)) call fatalError(Here, "Failed to get scattering reaction")

    ! Scatter particle
    if (reac % inCMFrame()) then
      collDat % A =  self % nuc % getMass()
      call self % scatterFromFixed(p, collDat, reac)
    else
      call self % scatterInLAB(p, collDat, reac)
    end if

    ! Apply weigth change
    p % w = p % w * reac % release(p % E)

  end subroutine inelastic

  !!
  !! Apply cutoffs
  !!
  subroutine cutoffs(self, p, tally, collDat, thisCycle, nextCycle)
    class(neutronCEstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(tallyAdmin), intent(inout)      :: tally
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle

    if (p % E < self % minE ) p % isDead = .true.

  end subroutine cutoffs

  !!
  !! Subroutine to perform scattering in LAB frame
  !! Returns mu -> cos of deflection angle in LAB frame
  !!
  subroutine scatterInLAB(self, p, collDat, reac)
    class(neutronCEstd), intent(inout)        :: self
    class(particle), intent(inout)            :: p
    type(collisionData), intent(inout)        :: collDat
    class(uncorrelatedReactionCE), intent(in) :: reac
    real(defReal)                             :: phi    ! Azimuthal scatter angle
    real(defReal)                             :: E_out, mu

    ! Sample scattering angles and post-collision energy
    call reac % sampleOut(mu, phi, E_out, p % E, p % pRNG)

    ! Update neutron state
    p % E = E_out
    call p % rotate(mu, phi)
    collDat % muL = mu

  end subroutine scatterInLAB

  !!
  !! Subroutine to perform scattering from stationary target.
  !! Returns mu -> cos of deflection angle in LAB frame
  !!
  subroutine scatterFromFixed(self, p, collDat, reac)
    class(neutronCEstd), intent(inout)         :: self
    class(particle), intent(inout)             :: p
    type(collisionData), intent(inout)         :: collDat
    class(uncorrelatedReactionCE), intent(in)  :: reac
    real(defReal)                              :: phi
    real(defReal)                              :: E_out
    real(defReal)                              :: E_outCM, mu
    integer(shortInt)                          :: MT

    ! Read data
    MT = collDat % MT

    ! Sample mu, phi and outgoing energy
    call reac % sampleOut(mu, phi, E_outCM, p % E, p % pRNG)

    ! Save incident energy
    E_out = p % E

    if (MT == N_N_elastic) then
      call asymptoticScatter(E_out, mu, collDat % A)
    else
      call asymptoticInelasticScatter(E_out, mu, E_outCM, collDat % A)
    end if

    ! Update particle state
    call p % rotate(mu, phi)
    p % E = E_out
    collDat % muL = mu

  end subroutine scatterFromFixed

  !!
  !! Subroutine to perform scattering from moving target
  !! Supports only elastic collisions
  !!
  subroutine scatterFromMoving(self, p, collDat, reac)
    class(neutronCEstd), intent(inout)         :: self
    class(particle), intent(inout)             :: p
    type(collisionData),intent(inout)          :: collDat
    class(uncorrelatedReactionCE), intent(in)  :: reac
    class(ceNeutronNuclide), pointer           :: ceNuc0K
    integer(shortInt)                          :: nucIdx
    real(defReal)                              :: A, kT, mu
    real(defReal),dimension(3)                 :: V_n           ! Neutron velocity (vector)
    real(defReal)                              :: U_n           ! Neutron speed (scalar)
    real(defReal),dimension(3)                 :: dir_pre       ! Pre-collision direction
    real(defReal),dimension(3)                 :: dir_post      ! Post-collicion direction
    real(defReal),dimension(3)                 :: V_t, V_cm     ! Target and CM velocity
    real(defReal)                              :: phi, dummy
    real(defReal)                              :: maj
    logical(defBool)                           :: inEnergyRange, hasDBRC
    character(100), parameter :: Here = 'ScatterFromMoving (neutronCEstd_class.f90)'

    ! Read collision data
    A      = collDat % A
    kT     = collDat % kT
    nucIdx = collDat % nucIdx

    ! Get neutron direction and velocity
    dir_pre = p % dirGlobal()
    V_n     = dir_pre * sqrt(p % E)

    ! Sample target velocity with constant XS or with DBRC
    ! Check energy range
    inEnergyRange = ((p % E <= self % DBRCeMax) .and. (self % DBRCeMin <= p % E))
    ! Check if DBRC is on for this target nuclide
    hasDBRC = self % nuc % hasDBRC()

    if (inEnergyRange .and. hasDBRC) then

      ! Retrieve 0K nuclide index from DBRC nuclide map
      nucIdx = self % xsData % mapDBRCnuc % get(nucIdx)

      ! Assign pointer for the 0K nuclide
      ceNuc0K => ceNeutronNuclide_CptrCast(self % xsData % getNuclide(nucIdx))
      if(.not.associated(ceNuc0K)) call fatalError(Here, 'Failed to retrieve CE Neutron Nuclide')

      ! Get elastic scattering 0K majorant
      maj = self % xsData % getScattMicroMajXS(p % E, kT, A, nucIdx)

      ! Use DBRC to sample target velocity
      V_t = targetVelocity_DBRCXS(ceNuc0K, p % E, dir_pre, A, kT, p % pRNG, maj)

    else
      ! Constant cross section approximation
      V_t = targetVelocity_constXS(p % E, dir_pre, A, kT, p % pRNG)

    end if

    ! Calculate Centre-of-Mass velocity
    V_cm = (V_n + V_t *A)/(A+1)

    ! Move Neutron velocity to CM frame, store speed and calculate new normalised direction
    V_n = V_n - V_cm
    U_n = norm2(V_n)
    V_n = V_n / U_n

    ! Sample mu and phi in CM frame
    call reac % sampleOut(mu, phi, dummy, p % E, p % pRNG)

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
