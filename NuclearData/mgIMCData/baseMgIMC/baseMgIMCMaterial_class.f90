module baseMgIMCMaterial_class

  use numPrecision
  use endfConstants
  use universalVariables
  use genericProcedures, only : fatalError, numToChar
  use RNG_class,         only : RNG
  use dictionary_class,  only : dictionary
  use poly_func

  ! Nuclear Data Interfaces
  use materialHandle_inter,    only : materialHandle
  use mgIMCMaterial_inter,     only : mgIMCMaterial, kill_super => kill
  use IMCXSPackages_class,     only : IMCMacroXSs
  use materialEquations

  use simulationTime_class,    only : timeStep

  implicit none
  private

  !!
  !! Public Pointer Cast
  !!
  public :: baseMgIMCMaterial_TptrCast
  public :: baseMgIMCMaterial_CptrCast

  ! Public data location parameters
  ! Use them if accessing data entries directly
  integer(shortInt), parameter, public :: TOTAL_XS      = 1
  integer(shortInt), parameter, public :: IESCATTER_XS  = 2
  integer(shortInt), parameter, public :: CAPTURE_XS    = 3
  integer(shortInt), parameter, public :: EMISSION_PROB = 4

  ! Calculation Type
  integer(shortInt), parameter, public :: IMC  = 1
  integer(shortInt), parameter, public :: ISMC = 2

  !!
  !! Material class for MG IMC and ISMC calculations
  !!
  !! Stores MG data in a table.
  !!
  !! Public members:
  !!   data -> Rank 2 array with all XSs data
  !!
  !! Interface:
  !!   materialHandle interface
  !!   mgIMCMaterial interface
  !!   init                 -> initialise Basic MG Material from dictionary and config keywords
  !!   nGroups              -> returns number of energy groups
  !!   updateMat            -> update material properties as required for IMC calculation  
  !!   getEmittedRad        -> returns the radiation to be emitted in current timestep
  !!   getFleck             -> returns current material Fleck factor
  !!   getEta               -> returns current value of eta (ISMC only)
  !!   getTemp              -> returns current material temperature
  !!   getMatEnergy         -> returns energy of material
  !!   setCalcType          -> set to IMC or ISMC
  !!   sampleEnergyGroup    -> return sampled energy group of an emitted photon
  !!   sampleTransformTime  -> return sampled time for transform of P_MATERIAL to P_PHOTON (ISMC)
  !!
  !! Note:
  !!   Order of "data" array is: data(XS_type, Group #)
  !!   Dictionary with data must contain following entries:
  !!     -> numberOfGroups
  !!
  type, public, extends(mgIMCMaterial) :: baseMgIMCMaterial
    real(defReal),dimension(:,:), allocatable :: data
    character(nameLen)                        :: name   ! Name for update equations (see materialEquations.f90)
    real(defReal)                             :: T      ! Temperature
    real(defReal)                             :: V      ! Volume
    real(defReal)                             :: fleck  ! Fleck factor
    real(defReal)                             :: alpha  ! User-defined parameter for fleck factor
    real(defReal)                             :: sigmaP ! Planck opacity
    real(defReal)                             :: matEnergy      ! Total energy stored in material
    real(defReal)                             :: prevMatEnergy  ! Energy prior to material update
    real(defReal)                             :: energyDens     ! Energy density = matEnergy/V
    real(defReal)                             :: eta            ! aT^4/energyDens, used for ISMC only
    integer(shortInt)                         :: calcType       ! IMC or ISMC
    real(defReal)                             :: sigmaFactor
    real(defReal)                             :: cvFactor

  contains
    ! Superclass procedures
    procedure :: kill
    procedure :: getMacroXSs_byG
    procedure :: getTotalXS

    ! Local procedures
    procedure :: init
    procedure :: nGroups
    procedure :: updateMat
    procedure :: getEmittedRad
    procedure :: getFleck
    procedure :: getEta
    procedure :: getTemp
    procedure :: getMatEnergy
    procedure :: setCalcType
    procedure :: sampleEnergyGroup
    procedure :: sampleTransformTime

    procedure, private :: tempFromEnergy
    procedure, private :: sigmaFromTemp
    procedure, private :: updateFleck

  end type baseMgIMCMaterial

contains

!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Standard procedures
!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(baseMgIMCMaterial), intent(inout) :: self

    ! Call superclass procedure
    call kill_super(self)

    ! Kill local content
    if(allocated(self % data)) deallocate(self % data)

  end subroutine kill

  !!
  !! Load Macroscopic XSs into the provided package for a given group index G
  !!
  !! See mgIMCMaterial documentation for more details
  !!
  subroutine getMacroXSs_byG(self, xss, G, rand)
    class(baseMgIMCMaterial), intent(in)     :: self
    type(IMCMacroXSs), intent(out)           :: xss
    integer(shortInt), intent(in)            :: G
    class(RNG), intent(inout)                :: rand
    character(100), parameter :: Here = ' getMacroXSs (baseMgIMCMaterial_class.f90)'

    ! Verify bounds
    if(G < 1 .or. self % nGroups() < G) then
      call fatalError(Here,'Invalid group number: '//numToChar(G)// &
                           ' Data has only: ' // numToChar(self % nGroups()))
    end if

    ! Get XSs
    xss % total            = self % data(TOTAL_XS, G)
    xss % elasticScatter   = ZERO
    xss % inelasticScatter = self % data(IESCATTER_XS, G)
    xss % capture          = self % data(CAPTURE_XS, G)

  end subroutine getMacroXSs_byG

  !!
  !! Return Total XSs for energy group G
  !!
  !! See mgIMCMaterial documentationfor details
  !!
  function getTotalXS(self, G, rand) result(xs)
    class(baseMgIMCMaterial), intent(in)     :: self
    integer(shortInt), intent(in)            :: G
    class(RNG), intent(inout)                :: rand
    real(defReal)                            :: xs
    character(100), parameter :: Here = ' getTotalXS (baseMgIMCMaterial_class.f90)'

    ! Verify bounds
    if(G < 1 .or. self % nGroups() < G) then
      call fatalError(Here,'Invalid group number: '//numToChar(G)// &
                           ' Data has only: ' // numToChar(self % nGroups()))
      xs = ZERO ! Avoid warning
    end if
    xs = self % data(TOTAL_XS, G)

  end function getTotalXS

  !!
  !! Initialise Base MG IMC Material fromdictionary
  !!
  !! Args:
  !!   dict       [in] -> Input dictionary with all required XSs
  !!
  !! Errors:
  !!   FatalError if data in dictionary is invalid (inconsistant # of groups;
  !!     -ve entries in P0 XSs)
  !!
  subroutine init(self, dict)
    class(baseMgIMCMaterial), intent(inout)     :: self
    class(dictionary),target, intent(in)        :: dict
    integer(shortInt)                           :: nG, N, i
    real(defReal)                               :: dT, tempT, tempU
    character(100), parameter :: Here = 'init (baseMgIMCMaterial_class.f90)'

    ! Read number of groups
    if (associated(mgEnergyGrid)) then
      nG = mgEnergyGrid % getSize()
      if (nG < 1) call fatalError(Here,'Number of groups is invalid' // numToChar(nG))
    else
      nG = 1
    end if

    ! Allocate space for data
    N = 4
    allocate(self % data(N, nG))

    ! Store alpha setting
    call dict % getOrDefault(self % alpha, 'alpha', ONE)

    ! Get name of equations for heat capacity and opacity calculations
    call dict % get(self % name, 'equations')

    ! Get optional multiplication factor for heat capacity and opacity
    call dict % getOrDefault(self % sigmaFactor, 'sigmaMultiple', ONE)
    call dict % getOrDefault(self % cvFactor, 'cvMultiple', ONE)

    ! Read initial temperature and volume
    call dict % get(self % T, 'T')
    call dict % get(self % V, 'V')

    ! Calculate initial opacities
    call self % sigmaFromTemp()

    ! Calculate initial energy density by integration
    dT = self % T / 1000
    tempT = dT/2
    tempU = 0
    do i=1, 1000
      tempU = tempU + dT * evaluateCv(self % name, tempT)
      tempT = tempT + dT
    end do
    self % energyDens = tempU * self % cvFactor
    self % matEnergy  = self % energyDens * self % V

  end subroutine init

  !!
  !! Return number of energy groups
  !!
  !! Args:
  !!   None
  !!
  !! Errors:
  !!   None
  !!
  pure function nGroups(self) result(nG)
    class(baseMgIMCMaterial), intent(in)     :: self
    integer(shortInt)                        :: nG

    if(allocated(self % data)) then
      nG = size(self % data,2)
    else
      nG = 0
    end if

  end function nGroups

  !!
  !! Cast materialHandle pointer to baseMgIMCMaterial type pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class materialHandle
  !!
  !! Result:
  !!   Null if source is not of baseMgIMCMaterial type
  !!   Target points to source if source is baseMgIMCMaterialtype
  !!
  pure function baseMgIMCMaterial_TptrCast(source) result(ptr)
    class(materialHandle), pointer, intent(in) :: source
    type(baseMgIMCMaterial), pointer           :: ptr

    select type(source)
      type is(baseMgIMCMaterial)
        ptr => source

      class default
        ptr => null()
    end select

  end function baseMgIMCMaterial_TptrCast

  !!
  !! Cast materialHandle pointer to baseMgIMCMaterial class pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class materialHandle
  !!
  !! Result:
  !!   Null if source is not of baseMgIMCMaterial class
  !!   Target points to source if source is baseMgIMCMaterial class
  !!
  pure function baseMgIMCMaterial_CptrCast(source) result(ptr)
    class(materialHandle), pointer, intent(in) :: source
    class(baseMgIMCMaterial), pointer          :: ptr

    select type(source)
      class is(baseMgIMCMaterial)
        ptr => source

      class default
        ptr => null()
    end select

  end function baseMgIMCMaterial_CptrCast

  !!
  !! Set the calculation type to be used
  !!
  !! Current options:
  !!   IMC
  !!   ISMC
  !!
  !! Errors:
  !!   Unrecognised option
  !!
  subroutine setCalcType(self, calcType)
    class(baseMgIMCMaterial), intent(inout) :: self
    integer(shortInt), intent(in)           :: calcType
    character(100), parameter               :: Here = 'setCalcType (baseMgIMCMaterial_class.f90)'

    if(calcType /= IMC .and. calcType /= ISMC) call fatalError(Here, 'Invalid calculation type')

    self % calcType = calcType

    ! Set initial fleck factor
    call self % updateFleck()

  end subroutine setCalcType


!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Material Updates
!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Update material properties at each time step
  !! First update energy using simple balance, then solve for temperature,
  !!  then update temperature-dependent properties
  !!
  !! Args:
  !!   tallyEnergy [in] -> Energy absorbed into material
  !!   printUpdate [in, optional] -> Bool, if true then will print updates to screen
  !!
  subroutine updateMat(self, tallyEnergy, printUpdate)
    class(baseMgIMCMaterial),intent(inout)  :: self
    real(defReal), intent(in)               :: tallyEnergy
    logical(defBool), intent(in), optional  :: printUpdate
    character(100), parameter               :: Here = "updateMat (baseMgIMCMaterial_class.f90)"

    ! TODO: Print updates if requested

    ! Save previous energy
    self % prevMatEnergy = self % matEnergy

    ! Update material internal energy
    if (self % calcType == IMC) then
      self % matEnergy = self % matEnergy - self % getEmittedRad() + tallyEnergy
    else
      self % matEnergy = tallyEnergy
    end if

    self % energyDens = self % matEnergy / self % V

    ! Return if no change
    if (abs(self % matEnergy - self % prevMatEnergy) < 0.00001*self % prevMatEnergy) return

    ! Confirm new energy density is valid
    if (self % energyDens <= ZERO) call fatalError(Here, 'Energy density is not positive')

    ! Update material temperature
    self % T = self % tempFromEnergy()

    ! Update sigma
    call self % sigmaFromTemp()

    ! Update fleck factor
    call self % updateFleck()

  end subroutine updateMat

  !!
  !! Calculate material temperature from internal energy by integration
  !!
  !! Step up (or down if energy is lower than in previous step) through temperature in steps of dT
  !! and increment energy density by dT*cv(T) until target energy density is reached
  !!
  function tempFromEnergy(self) result(T)
    class(baseMgIMCMaterial), intent(inout) :: self
    real(defReal)                           :: T, dT, tempT, U, tempU, tol, error, increase
    integer(shortInt)                       :: i
    character(100), parameter               :: Here = 'tempFromEnergy (mgIMCMaterial_class.f90)'

    ! Parameters affecting accuracy
    dT  = self % T / 1000     ! Initial temperature step size to take, reduced after overshoot
    tol = self % T / 1000000  ! Continue stepping until within tolerance

    ! Starting temperature and energy density
    T = self % T
    U = self % prevMatEnergy / self % V

    ! If starting energy density is higher than target, flip dT to be negative
    if (U > self % energyDens) dT = -dT

    i = 0

    integrate:do

      ! Protect against infinite loop
      i = i+1
      if (i > 100000) then
        print *, 'Energy density: ', self % energyDens
        call fatalError(Here, "100,000 iterations without convergence, maybe NaN energy density?")
      end if
      ! Increase step size to avoid lack of convergence due to very small starting temperature
      if (mod(i,1000)==0) dT = 10*dT

      ! Increment temperature and increment the corresponding energy density
      tempT = T + dT/2
      increase = dT * evaluateCv(self % name, tempT) * self % cvFactor
      ! Protect against division by 0 or other numerical errors
      if (increase /= increase .or. increase > INF) increase = ZERO
      tempU = U + increase

      error = self % energyDens - tempU

      if (abs(error) < tol) then ! Finished
        T = T + dT
        exit integrate
      end if

      if (error*dT < 0) then ! If error and dT have different signs then we have overshot
        ! Decrease dT and try again
        dT = dT/2
        cycle integrate
      end if

      ! Update temperature and energy and return to start
      T = T + dT
      U = tempU

    end do integrate

  end function tempFromEnergy

  !!
  !! Calculate opacities from current temp
  !!
  subroutine sigmaFromTemp(self)
    class(baseMgIMCMaterial), intent(inout) :: self
    integer(shortInt)                       :: i, j, stepsPerGroup
    real(defReal)                           :: sigmaP, E, EStep, increase, upper, lower
    character(100), parameter :: Here = 'sigmaFromTemp (baseMgIMCMaterial_class.f90)'

    ! Evaluate opacities for grey case
    if (self % nGroups() == 1) then
      self % data(CAPTURE_XS,1)   = evaluateSigma(self % name, self % T, ONE) * self % sigmaFactor
      self % data(IESCATTER_XS,1) = ZERO
      self % data(TOTAL_XS,1)     = self % data(CAPTURE_XS,1) + self % data(IESCATTER_XS,1)
      ! Planck opacity equal to absorption opacity for single frequency
      self % sigmaP = self % data(CAPTURE_XS, 1)
      return
    end if

    ! Evaluate opacities for frequency-dependent case
    do i = 1, self % nGroups()
      ! Take geometric mean of upper and lower group boundaries
      upper = evaluateSigma(self % name, self % T, mgEnergyGrid % bin(i))
      lower = evaluateSigma(self % name, self % T, mgEnergyGrid % bin(i+1))
      self % data(CAPTURE_XS,i) = sqrt(upper*lower)*self % sigmaFactor !min(1e5_defReal, sqrt(o1 * o2))

      upper = upper * normPlanckSpectrum(mgEnergyGrid % bin(i), self % T)
      lower = lower * normPlanckSpectrum(mgEnergyGrid % bin(i+1), self % T)
      self % data(EMISSION_PROB,i) = sqrt(upper*lower)
      self % data(EMISSION_PROB,i) = self % data(EMISSION_PROB,i)*(mgEnergyGrid % bin(i) - mgEnergyGrid % bin(i+1))

    end do
    self % data(IESCATTER_XS,:) = ZERO
    self % data(TOTAL_XS,:)     = self % data(CAPTURE_XS,:) + self % data(IESCATTER_XS,:)
    self % data(EMISSION_PROB,:) = self % data(EMISSION_PROB,:)/sum(self % data(EMISSION_PROB,:))

    ! Evaluate Planck opacity via a much finer numerical integration
    sigmaP = ZERO
    do i = 1, self % nGroups()
      ! 100 integration steps per energy group seems to give very good accuracy, might be worth playing around with
      stepsPerGroup = 100
      EStep = (mgEnergyGrid % bin(i) - mgEnergyGrid % bin(i+1)) / stepsPerGroup
      E = mgEnergyGrid % bin(i) - 0.5*EStep
      do j = 1, stepsPerGroup
        increase = normPlanckSpectrum(E, self % T)*evaluateSigma(self % name, self % T, E)
        sigmaP = sigmaP + EStep * increase
        E = E - EStep
      end do
    end do
    self % sigmaP = sigmaP * self % sigmaFactor

  end subroutine sigmaFromTemp

  !!
  !! Update fleck factor
  !!
  subroutine updateFleck(self)
    class(baseMgIMCMaterial), intent(inout) :: self
    real(defReal)                           :: beta, zeta
    character(100), parameter               :: Here = 'updateFleck (baseMgIMCMaterial_class.f90)'

    ! Calculate beta, ratio of radiation and material heat capacities
    beta = 4 * radiationConstant * self % T**3 / (evaluateCv(self % name, self % T)*self % cvFactor)

    ! Use time step size to calculate fleck factor
    select case(self % calcType)

      case(IMC)
        self % fleck = 1/(1+self % sigmaP*lightSpeed*beta*timeStep()*self % alpha)

      case(ISMC)
        self % eta = radiationConstant * self % T**4 / self % energyDens
        zeta = beta - self % eta
        self % fleck = 1 / (1 + zeta*self % sigmaP*lightSpeed*timeStep())

      case default
        call fatalError(Here, 'Unrecognised calculation type')

    end select

  end subroutine updateFleck


!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Obtain material info
!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Return the energy to be emitted during time step, E_r
  !!
  function getEmittedRad(self) result(emittedRad)
    class(baseMgIMCMaterial), intent(inout) :: self
    real(defReal)                           :: u_r, emittedRad

    u_r = radiationConstant * (self % T)**4

    emittedRad = lightSpeed * timeStep() * self % sigmaP * self % fleck * u_r * self % V

  end function getEmittedRad

  !!
  !! Return the fleck factor of the material
  !!
  function getFleck(self) result(fleck)
    class(baseMgIMCMaterial),intent(in) :: self
    real(defReal)                       :: fleck

    fleck = self % fleck

  end function getFleck

  !!
  !! Return eta = aT**4/U_m
  !!
  !! Currently only used in transportOperatorIMC_class.f90 for ISMC calculations
  !!
  function getEta(self) result(eta)
    class(baseMgIMCMaterial),intent(in) :: self
    real(defReal)                       :: eta

    eta = self % eta

  end function getEta

  !!
  !! Get temperature of material
  !!
  function getTemp(self) result(T)
    class(baseMgIMCMaterial), intent(inout) :: self
    real(defReal)                           :: T

    T = self % T

  end function getTemp

  !!
  !! Return energy per unit volume of material
  !!
  function getMatEnergy(self) result(energy)
    class(baseMgIMCMaterial), intent(inout) :: self
    real(defReal)                           :: energy

    energy = self % matEnergy

  end function getMatEnergy


!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Sample emission properties
!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Sample energy group for a particle emitted from material (and after effective scattering)
  !!
  function sampleEnergyGroup(self, rand) result(G)
    class(baseMgIMCMaterial), intent(inout) :: self
    class(RNG), intent(inout)               :: rand
    integer(shortInt)                       :: G, i
    real(defReal)                           :: random, cumProb

    if (self % nGroups() == 1) then
      G = 1
      return
    end if

    random  = rand % get()
    cumProb = ZERO

    ! Choose based on emission probability of each group
    do i = 1, self % nGroups()
      cumProb = cumProb + self % data(EMISSION_PROB,i)
      if (random < cumProb) then
        G = i
        return
      end if
    end do

  end function sampleEnergyGroup   

  !!
  !! Sample the time taken for a material particle to transform into a photon
  !! Used for ISMC only
  !!
  function sampleTransformTime(self, rand) result(t)
    class(baseMgIMCMaterial), intent(inout) :: self
    class(RNG), intent(inout)               :: rand
    real(defReal)                           :: t

    t = -log(rand % get()) / (self % sigmaP * self % fleck * self % eta * lightSpeed)

  end function sampleTransformTime

end module baseMgIMCMaterial_class
