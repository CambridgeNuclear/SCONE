module neutronCEimp_class

  use numPrecision
  use endfConstants
  use universalVariables,            only : nameUFS, nameWW
  use genericProcedures,             only : fatalError, rotateVector, numToChar
  use dictionary_class,              only : dictionary
  use RNG_class,                     only : RNG

  ! Particle types
  use particle_class,                only : particle, particleState, printType, P_NEUTRON
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

  ! Geometry and fields
  use geometryReg_mod,                only : gr_fieldIdx => fieldIdx, gr_fieldPtr => fieldPtr
  use uniFissSitesField_class,        only : uniFissSitesField, uniFissSitesField_TptrCast
  use weightWindowsField_class,       only : weightWindowsField, weightWindowsField_TptrCast

  ! Cross-Section Packages
  use neutronXsPackages_class,       only : neutronMicroXSs

  ! Scattering procedures
  use scatteringKernels_func, only : asymptoticScatter, targetVelocity_constXS, &
                                     asymptoticInelasticScatter
  implicit none
  private

  !!
  !! Standard scalar collision processor for CE neutrons
  !!   -> Preforms implicit or analog fission site generation
  !!   -> Preforms implicit or analog capture
  !!   -> Treats fission as capture (only implicit generation of 2nd-ary neutrons)
  !!   -> Does not create secondary non-neutron projectiles
  !!
  !! Settings:
  !!  minE    -> minimum energy cut-off [MeV] (default = 1.0E-11)
  !!  maxE    -> maximum energy. Higher energies are set to maximum (not re-rolled) [MeV]
  !!             (default = 20.0)
  !!  minWgt  -> minimum particle weight for rouletting (optional)
  !!  maxWgt  -> maximum particle weight for splitting (optional)
  !!  avgWgt  -> weight of a particle on surviving rouletting (optional)
  !!  impAbs  -> is implicit capture performed? (off by default)
  !!  impGen  -> are fission sites generated implicitly? (on by default)
  !!  UFS     -> uniform fission sites variance reduction
  !!  thresh_E -> Energy threshold for explicit treatment of target nuclide movement [-].
  !!              Target movment is sampled if neutron energy E < kT * thresh_E where
  !!              kT is target material temperature in [MeV]. (default = 400.0)
  !!  thresh_A -> Mass threshold for explicit tratment of target nuclide movement [Mn].
  !!              Target movment is sampled if target mass A < thresh_A. (default = 1.0)
  !!  splitting -> splits particles above certain weight (on by default)
  !!  roulette  -> roulettes particles below certain weight (off by defautl)
  !!  weightWindows -> uses a weight windows field (off by default)
  !!
  !! Sample dictionary input:
  !!   collProcName {
  !!   type            neutronCEimp;
  !!   #minEnergy      <real>;#
  !!   #maxEnergy      <real>;#
  !!   #energyTreshold <real>;#
  !!   #massTreshold   <real>;#
  !!   #splitting      <logical>;#
  !!   #roulette       <logical>;#
  !!   #minWgt         <real>;#
  !!   #maxWgt         <real>;#
  !!   #avgWgt         <real>;#
  !!   #impAbs         <logical>;#
  !!   #impGen         <logical>;#
  !!   #UFS            <logical>;#
  !!   #weightWindows  <logical>;#
  !!   }
  !!
  type, public, extends(collisionProcessor) :: neutronCEimp
    private
    !! Nuclear Data block pointer -> public so it can be used by subclasses (protected member)
    class(ceNeutronDatabase), pointer, public :: xsData => null()
    class(ceNeutronMaterial), pointer, public :: mat    => null()
    class(ceNeutronNuclide),  pointer, public :: nuc    => null()
    class(uniFissSitesField), pointer :: ufsField => null()

    !! Settings - private
    real(defReal) :: minE
    real(defReal) :: maxE
    real(defReal) :: minWgt
    real(defReal) :: maxWgt
    real(defReal) :: avWgt
    real(defReal) :: thresh_E
    real(defReal) :: thresh_A
    ! Variance reduction options
    logical(defBool) :: weightWindows
    logical(defBool) :: splitting
    logical(defBool) :: roulette
    logical(defBool) :: implicitAbsorption ! Prevents particles dying through capture
    logical(defBool) :: implicitSites ! Generates fission sites on every fissile collision
    logical(defBool) :: uniFissSites

    type(weightWindowsField), pointer :: weightWindowsMap

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

    ! Variance reduction procedures
    procedure, private :: split
    procedure, private :: russianRoulette
  end type neutronCEimp

contains

  !!
  !! Initialise from dictionary
  !!
  subroutine init(self, dict)
    class(neutronCEimp), intent(inout) :: self
    class(dictionary), intent(in)      :: dict
    integer(shortInt)                  :: idx
    character(100), parameter :: Here = 'init (neutronCEimp_class.f90)'

    ! Call superclass
    call init_super(self, dict)

    ! Read settings for neutronCEimp
    ! Maximum and minimum energy
    call dict % getOrDefault(self % minE,'minEnergy',1.0E-11_defReal)
    call dict % getOrDefault(self % maxE,'maxEnergy',20.0_defReal)

    ! Thermal scattering kernel thresholds
    call dict % getOrDefault(self % thresh_E, 'energyThreshold', 400.0_defReal)
    call dict % getOrDefault(self % thresh_A, 'massThreshold', 1.0_defReal)

    ! Obtain settings for variance reduction
    call dict % getOrDefault(self % weightWindows,'weightWindows', .false.)
    call dict % getOrDefault(self % splitting,'split', .false.)
    call dict % getOrDefault(self % roulette,'roulette', .false.)
    call dict % getOrDefault(self % minWgt,'minWgt',0.25_defReal)
    call dict % getOrDefault(self % maxWgt,'maxWgt',1.25_defReal)
    call dict % getOrDefault(self % avWgt,'avWgt',0.5_defReal)
    call dict % getOrDefault(self % implicitAbsorption,'impAbs', .false.)
    call dict % getOrDefault(self % implicitSites,'impGen', .true.)
    call dict % getOrDefault(self % uniFissSites,'UFS', .false.)

    ! Verify settings
    if( self % minE < ZERO ) call fatalError(Here,'-ve minEnergy')
    if( self % maxE < ZERO ) call fatalError(Here,'-ve maxEnergy')
    if( self % minE >= self % maxE) call fatalError(Here,'minEnergy >= maxEnergy')
    if( self % thresh_E < 0) call fatalError(Here,' -ve energyThreshold')
    if( self % thresh_A < 0) call fatalError(Here,' -ve massThreshold')

    if (self % splitting) then
      if (self % maxWgt < 2 * self % minWgt) call fatalError(Here,&
              'Upper weight bound must be at least twice the lower weight bound')
    end if

    if (self % implicitAbsorption) then
      if (.not.self % roulette) call fatalError(Here,&
         'Must use Russian roulette when using implicit absorption')
      if (.not.self % implicitSites) call fatalError(Here,&
         'Must generate fission sites implicitly when using implicit absorption')
    end if

    ! Sets up the uniform fission sites field
    if (self % uniFissSites) then
      idx = gr_fieldIdx(nameUFS)
      self % ufsField => uniFissSitesField_TptrCast(gr_fieldPtr(idx))
    end if

    ! Sets up the weight windows field
    if (self % weightWindows) then
      idx = gr_fieldIdx(nameWW)
      self % weightWindowsMap => weightWindowsField_TptrCast(gr_fieldPtr(idx))
    end if

  end subroutine init

  !!
  !! Samples collision without any implicit treatment
  !!
  subroutine sampleCollision(self, p, collDat, thisCycle, nextCycle)
    class(neutronCEimp), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle
    type(neutronMicroXSs)                :: microXSs
    real(defReal)                        :: r
    character(100),parameter :: Here = 'sampleCollision (neutronCEimp_class.f90)'

    ! Verify that particle is CE neutron
    if(p % isMG .or. p % type /= P_NEUTRON) then
      call fatalError(Here, 'Supports only CE Neutron. Was given MG '//printType(p % type))
    end if

    ! Verify and load nuclear data pointer
    self % xsData => ndReg_getNeutronCE()
    if(.not.associated(self % xsData)) call fatalError(Here, 'There is no active Neutron CE data!')

    ! Verify and load material pointer
    self % mat => ceNeutronMaterial_CptrCast( self % xsData % getMaterial( p % matIdx()))
    if(.not.associated(self % mat)) call fatalError(Here, 'Material is not ceNeutronMaterial')

    ! Select collision nuclide
    collDat % nucIdx = self % mat % sampleNuclide(p % E, p % pRNG)

    self % nuc => ceNeutronNuclide_CptrCast(self % xsData % getNuclide(collDat % nucIdx))
    if(.not.associated(self % mat)) call fatalError(Here, 'Failed to retive CE Neutron Nuclide')

    ! Select Main reaction channel
    call self % nuc % getMicroXSs(microXss, p % E, p % pRNG)
    r = p % pRNG % get()
    collDat % MT = microXss % invert(r)

  end subroutine sampleCollision

  !!
  !! Perform implicit treatment
  !!
  subroutine implicit(self, p, collDat, thisCycle, nextCycle)
    class(neutronCEimp), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle
    type(fissionCE), pointer             :: fission
    type(neutronMicroXSs)                :: microXSs
    type(particleState)                  :: pTemp
    real(defReal),dimension(3)           :: r, dir, val
    integer(shortInt)                    :: n, i
    real(defReal)                        :: wgt, rand1, E_out, mu, phi
    real(defReal)                        :: sig_nufiss, sig_tot, k_eff, &
                                            sig_scatter, totalElastic
    logical(defBool)                     :: fiss_and_implicit
    character(100),parameter             :: Here = 'implicit (neutronCEimp_class.f90)'

    ! Generate fission sites if nuclide is fissile
    fiss_and_implicit = self % nuc % isFissile() .and. self % implicitSites
    if (fiss_and_implicit) then
      ! Obtain required data
      wgt   = p % w                ! Current weight
      k_eff = p % k_eff            ! k_eff for normalisation
      rand1 = p % pRNG % get()     ! Random number to sample sites

      call self % nuc % getMicroXSs(microXSs, p % E, p % pRNG)
      sig_nufiss = microXSs % nuFission
      sig_tot    = microXSs % total

      ! Sample number of fission sites generated
      ! Support -ve weight particles
      if (self % uniFissSites) then
        val = self % ufsField % at(p)
        n = int(abs( (wgt * sig_nufiss) / (sig_tot * k_eff))*val(1)/val(2) + rand1, shortInt)
        wgt =  val(2)/val(1)
      else
        n = int(abs( (wgt * sig_nufiss) / (sig_tot * k_eff)) + rand1, shortInt)
        wgt =  sign(ONE, wgt)
      end if

      ! Shortcut particle generation if no particles were sampled
      if (n < 1) return

      ! Get fission Reaction
      fission => fissionCE_TptrCast(self % xsData % getReaction(N_FISSION, collDat % nucIdx))
      if(.not.associated(fission)) call fatalError(Here, "Failed to get fissionCE")

      ! Store new sites in the next cycle dungeon
      r   = p % rGlobal()

      do i=1,n
        call fission % sampleOut(mu, phi, E_out, p % E, p % pRNG)
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

        call nextCycle % detain(pTemp)
        if (self % uniFissSites) call self % ufsField % storeFS(pTemp)

      end do
    end if

    ! Perform implicit absorption
    if (self % implicitAbsorption) then
      if(.not.fiss_and_implicit) then
        call self % nuc % getMicroXSs(microXSs, p % E, p % pRNG)
      end if
      sig_scatter  = microXSs % elasticScatter + microXSs % inelasticScatter
      sig_tot      = microXSs % total
      p % w        = p % w * sig_scatter/sig_tot
      ! Sample between elastic and inelastic
      totalElastic = microXSs % elasticScatter + microXSs % inelasticScatter
      if (p % pRNG % get() < microXSs % elasticScatter/totalElastic) then
        collDat % MT = N_N_elastic
      else
        collDat % MT = N_N_inelastic
      end if
    end if

  end subroutine implicit

  !!
  !! Process capture reaction
  !!
  subroutine capture(self, p, collDat, thisCycle, nextCycle)
    class(neutronCEimp), intent(inout)   :: self
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
    class(neutronCEimp), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle
    type(neutronMicroXSs)                :: microXSs
    type(fissionCE), pointer             :: fiss
    type(particleState)                  :: pTemp
    real(defReal),dimension(3)           :: r, dir, val
    integer(shortInt)                    :: n, i
    real(defReal)                        :: wgt, rand1, E_out, mu, phi
    real(defReal)                        :: sig_nufiss, sig_fiss, k_eff
    character(100),parameter             :: Here = 'fission (neutronCEimp_class.f90)'

    if (.not.self % implicitSites) then
      ! Obtain required data
      wgt   = p % w                ! Current weight
      k_eff = p % k_eff            ! k_eff for normalisation
      rand1 = p % pRNG % get()     ! Random number to sample sites

      call self % nuc % getMicroXSs(microXSs, p % E, p % pRNG)
      sig_nufiss = microXSs % nuFission
      sig_fiss   = microXSs % fission

      ! Sample number of fission sites generated
      ! Support -ve weight particles
      ! Note change of denominator (sig_fiss) wrt implicit generation
      if (self % uniFissSites) then
        val = self % ufsField % at(p)
        n = int(abs( (wgt * sig_nufiss) / (sig_fiss * k_eff))*val(1)/val(2) + rand1, shortInt)
        wgt =  val(2)/val(1)
      else
        n = int(abs( (wgt * sig_nufiss) / (sig_fiss * k_eff)) + rand1, shortInt)
        wgt =  sign(ONE, wgt)
      end if

      ! Shortcut particle generation if no particles were sampled
      if (n < 1) return

      ! Get fission Reaction
      fiss => fissionCE_TptrCast(self % xsData % getReaction(N_FISSION, collDat % nucIdx))
      if(.not.associated(fiss)) call fatalError(Here, "Failed to get fissionCE")

      ! Store new sites in the next cycle dungeon
      r   = p % rGlobal()

      do i=1,n
        call fiss % sampleOut(mu, phi, E_out, p % E, p % pRNG)
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

        call nextCycle % detain(pTemp)
        if (self % uniFissSites) call self % ufsField % storeFS(pTemp)

      end do
    end if

    p % isDead =.true.

  end subroutine fission

  !!
  !! Process elastic scattering
  !!
  !! All CE elastic scattering happens in the CM frame
  !!
  subroutine elastic(self, p, collDat, thisCycle, nextCycle)
    class(neutronCEimp), intent(inout)     :: self
    class(particle), intent(inout)         :: p
    type(collisionData), intent(inout)     :: collDat
    class(particleDungeon),intent(inout)   :: thisCycle
    class(particleDungeon),intent(inout)   :: nextCycle
    class(uncorrelatedReactionCE), pointer :: reac
    logical(defBool)                       :: isFixed
    character(100),parameter :: Here = 'elastic (neutronCEimp_class.f90)'

    ! Get reaction
    reac => uncorrelatedReactionCE_CptrCast( self % xsData % getReaction(collDat % MT, collDat % nucIdx))
    if(.not.associated(reac)) call fatalError(Here,'Failed to get elastic neutron scatter')

    ! Scatter particle
    collDat % A =  self % nuc % getMass()
    collDat % kT = self % nuc % getkT()

    isFixed = (p % E > collDat % kT * self % thresh_E) .and. (collDat % A > self % thresh_A)

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
  subroutine inelastic(self, p, collDat, thisCycle, nextCycle)
    class(neutronCEimp), intent(inout)     :: self
    class(particle), intent(inout)         :: p
    type(collisionData), intent(inout)     :: collDat
    class(particleDungeon),intent(inout)   :: thisCycle
    class(particleDungeon),intent(inout)   :: nextCycle
    class(uncorrelatedReactionCE), pointer :: reac
    character(100),parameter  :: Here =' inelastic (neutronCEimp_class.f90)'

    ! Invert inelastic scattering and Get reaction
    collDat % MT = self % nuc % invertInelastic(p % E, p % pRNG)
    reac => uncorrelatedReactionCE_CptrCast( self % xsData % getReaction(collDat % MT, collDat % nucIdx))
    if(.not.associated(reac)) call fatalError(Here, "Failed to get scattering reaction")

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
  subroutine cutoffs(self, p, collDat, thisCycle, nextCycle)
    class(neutronCEimp), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle
    real(defReal), dimension(3)          :: val
    real(defReal)                        :: minWgt, maxWgt, avWgt

    if (p % E < self % minE) then
      p % isDead = .true.

    ! Weight Windows treatment
    elseif (self % weightWindows) then
      val = self % weightWindowsMap % at(p)
      minWgt = val(1)
      maxWgt = val(2)
      avWgt  = val(3)

      ! If a particle is outside the WW map and all the weight limits
      ! are zero nothing happens. NOTE: this holds for positive weights only
      if ((p % w > maxWgt) .and. (maxWgt /= ZERO)) then
        call self % split(p, thisCycle, maxWgt)
      elseif (p % w < minWgt) then
        call self % russianRoulette(p, avWgt)
      end if

    ! Splitting with fixed threshold
    elseif ((self % splitting) .and. (p % w > self % maxWgt)) then
      call self % split(p, thisCycle, self % maxWgt)

    ! Roulette with fixed threshold and survival weight
    elseif ((self % roulette) .and. (p % w < self % minWgt)) then
      call self % russianRoulette(p, self % avWgt)

    end if

  end subroutine cutoffs

  !!
  !! Perform Russian roulette on a particle
  !!
  subroutine russianRoulette(self, p, avWgt)
    class(neutronCEimp), intent(inout) :: self
    class(particle), intent(inout)     :: p
    real(defReal), intent(in)          :: avWgt

    if (p % pRNG % get() < (ONE - p % w/avWgt)) then
      p % isDead = .true.
    else
      p % w = avWgt
    end if

  end subroutine russianRoulette

  !!
  !! Split particle which has too large a weight
  !!
  subroutine split(self, p, thisCycle, maxWgt)
    class(neutronCEimp), intent(inout)    :: self
    class(particle), intent(inout)        :: p
    class(particleDungeon), intent(inout) :: thisCycle
    real(defReal), intent(in)             :: maxWgt
    integer(shortInt)                     :: mult, i

    ! This value must be at least 2
    mult = ceiling(p % w/maxWgt)

    ! Decrease weight
    p % w = p % w/mult

    ! Add split particle's to the dungeon
    do i = 1,mult-1
      call thisCycle % detain(p)
    end do

  end subroutine split

  !!
  !! Subroutine to perform scattering in LAB frame
  !! Returns mu -> cos of deflection angle in LAB frame
  !!
  subroutine scatterInLAB(self, p, collDat, reac)
    class(neutronCEimp), intent(inout)        :: self
    class(particle), intent(inout)            :: p
    type(collisionData), intent(inout)        :: collDat
    class(uncorrelatedReactionCE), intent(in) :: reac
    real(defReal)                             :: phi    ! Azimuthal scatter angle
    real(defReal)                             :: E_out, mu
    integer(shortInt)                         :: MT, nucIdx

    ! Read data
    MT = collDat % MT
    nucIdx = collDat % nucIdx

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
    class(neutronCEimp), intent(inout)         :: self
    class(particle), intent(inout)             :: p
    type(collisionData), intent(inout)         :: collDat
    class(uncorrelatedReactionCE), intent(in)  :: reac
    real(defReal)                              :: phi
    real(defReal)                              :: E_out
    real(defReal)                              :: E_outCM, mu
    integer(shortInt)                          :: MT, nucIdx

    ! Read data
    MT     = collDat % MT
    nucIdx = collDat % nucIdx

    ! Sample mu , phi and outgoing energy
    call reac % sampleOut(mu, phi, E_outCM, p % E, p % pRNG)

    ! Save incident energy
    E_out = p % E

    if( MT == N_N_elastic) then
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
    class(neutronCEimp), intent(inout)         :: self
    class(particle), intent(inout)             :: p
    type(collisionData),intent(inout)          :: collDat
    class(uncorrelatedReactionCE), intent(in)  :: reac
    integer(shortInt)                          :: MT, nucIdx
    real(defReal)                              :: A, kT, mu
    real(defReal),dimension(3)                 :: V_n           ! Neutron velocity (vector)
    real(defReal)                              :: U_n           ! Neutron speed (scalar)
    real(defReal),dimension(3)                 :: dir_pre       ! Pre-collision direction
    real(defReal),dimension(3)                 :: dir_post      ! Post-collicion direction
    real(defReal),dimension(3)                 :: V_t, V_cm     ! Target and CM velocity
    real(defReal)                              :: phi, dummy

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


end module neutronCEimp_class
