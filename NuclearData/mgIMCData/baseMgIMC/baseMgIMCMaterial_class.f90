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
  integer(shortInt), parameter, public :: PLANCK_XS     = 4

  !!
  !! Basic type of MG material data
  !!
  !! Stores MG data in a table.
  !! All other scattering reactions are lumped into single multiplicative scattering,
  !! which is stored as INELASTIC scatering in macroXSs package! After all it is inelastic in
  !! the sense that outgoing group can change. Diffrent types of multiplicative scattering can be
  !! build. See doc of "init" procedure for details.
  !!
  !! Public members:
  !!   data -> Rank 2 array with all XSs data
  !!
  !! Interface:
  !!   materialHandle interface
  !!   mgIMCMaterial interface
  !!   init -> initialise Basic MG Material from dictionary and config keyword
  !!   nGroups -> returns number of energy groups
  !!   updateMat -> update material properties as required for IMC calculation
  !!   getEmittedRad -> returns the radiation to be emitted in current timestep
  !!   getFleck -> returns current material Fleck factor
  !!   getTemp -> returns current material temperature
  !!
  !! Note:
  !!   Order of "data" array is: data(XS_type, Group #)
  !!   Dictionary with data must contain following entries:
  !!     -> numberOfGroups
  !!
  type, public, extends(mgIMCMaterial) :: baseMgIMCMaterial
    real(defReal),dimension(:,:), allocatable :: data
    real(defReal),dimension(:), allocatable   :: cv
    real(defReal),dimension(:), allocatable   :: updateEqn
    real(defReal),dimension(:), allocatable   :: absEqn
    real(defReal),dimension(:), allocatable   :: scattEqn
    real(defReal),dimension(:), allocatable   :: planckEqn
    real(defReal)                             :: T
    real(defReal)                             :: V
    real(defReal)                             :: fleck
    real(defReal)                             :: alpha
    real(defReal)                             :: deltaT
    real(defReal)                             :: sigmaP
    real(defReal)                             :: matEnergy
    real(defReal)                             :: eta
    integer(shortInt)                         :: calcType

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
    procedure :: getEnergyDens
    procedure :: setType
    procedure :: setTimeStep

    procedure, private :: updateMatIMC
    procedure, private :: updateMatISMC
    procedure, private :: tempFromEnergy
    procedure, private :: sigmaFromTemp

  end type baseMgIMCMaterial

contains

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
    xss % planck           = self % data(PLANCK_XS, G)

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
    integer(shortInt)                           :: nG, N
    real(defReal), dimension(:), allocatable    :: temp
    character(100), parameter :: Here = 'init (baseMgIMCMaterial_class.f90)'

    ! Read number of groups
    call dict % get(nG, 'numberOfGroups')
    if(nG < 1) call fatalError(Here,'Number of groups is invalid' // numToChar(nG))

    ! Allocate space for data
    N = 4
    allocate(self % data(N, nG))

    ! Store alpha setting
    call dict % getOrDefault(self % alpha, 'alpha', ONE)

    ! Read opacity equations
    call dict % get(temp, 'sigmaA')
    self % absEqn = temp
    call dict % get(temp, 'sigmaS')
    self % scattEqn = temp

    ! Build planck opacity equation
    ! For grey case, sigmaP = sigmaA. Will become more complicated for frequency-dependent case
    self % planckEqn = self % absEqn

    ! Read heat capacity equation
    call dict % get(temp, 'cv')
    self % cv = temp

    ! Build update equation
    call poly_integrate(temp)
    self % updateEqn = temp

    ! Read initial temperature and volume
    call dict % get(self % T, 'T')
    call dict % get(self % V, 'V')

    ! Calculate initial opacities and energy
    call self % sigmaFromTemp()
    self % matEnergy = poly_eval(self % updateEqn, self % T) * self % V

    ! Set calculation type (will support ISMC in the future)
    self % calcType = IMC

  end subroutine init


  !!
  !! Provide material with time step size
  !!
  !! Args:
  !!   dt [in] -> time step size [s]
  !!
  !! Errors:
  !!   fatalError if calculation type is invalid (valid options are IMC or ISMC)
  !!
  subroutine setTimeStep(self, dt)
    class(baseMgIMCMaterial), intent(inout) :: self
    real(defReal), intent(in)               :: dt
    real(defReal)                           :: beta, zeta
    character(100), parameter               :: Here = 'setTimeStep (baseMgIMCMaterial_class.f90)'

    self % deltaT = dt

    ! Use time step size to calculate fleck factor
    if(self % calcType == IMC) then
      self % fleck = 1/(1+1*self % sigmaP*lightSpeed*self % deltaT*self % alpha)

    else if(self % calcType == ISMC) then
      beta = 4*radiationConstant * self % T**3 / poly_eval(self % cv, self % T)
      self % eta  =   radiationConstant * self % T**4 / self % matEnergy
      zeta = beta - self % eta
      self % fleck = 1 / (1 + zeta*self % sigmaP*lightSpeed*self % deltaT)

    else
      call fatalError(Here, 'Calculation type invalid or not set')
    end if


  end subroutine setTimeStep

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

    select case (self % calcType)

      case(IMC)
        call self % updateMatIMC(tallyEnergy, printUpdate)

      case(ISMC)
        call self % updateMatISMC(tallyEnergy, printUpdate)

      case default
        call fatalError(Here, "Invalid calculation type")

    end select

  end subroutine updateMat

  !!
  !! Material update for IMC calculation
  !!
  subroutine updateMatIMC(self, tallyEnergy, printUpdate)
    class(baseMgIMCMaterial), intent(inout) :: self
    real(defReal), intent(in)               :: tallyEnergy
    logical(defBool), intent(in), optional  :: printUpdate
    character(100), parameter               :: Here = "updateMatIMC (baseMgIMCMaterial_class.f90)"

    ! Print current properties
    if (present(printUpdate)) then
      if (printUpdate .eqv. .True.) then
        print *, "  T_old =                         ", self % T
        print *, "  matEnergy at start of timestep =", self % matEnergy
        print *, "  emittedRad =                    ", self % getEmittedRad()
        print *, "  tallyEnergy =                   ", tallyEnergy
      end if
    end if

    ! Update material internal energy
    self % matEnergy = self % matEnergy - self % getEmittedRad() + tallyEnergy

    ! Update material temperature
    self % T = self % tempFromEnergy()

    ! Update sigma
    call self % sigmaFromTemp()

    if( self % T < 0 ) then
     call fatalError(Here, "Temperature is negative")
    end if

    self % fleck = 1/(1+1*self % sigmaP*lightSpeed*self % deltaT*self % alpha)

    ! Print updated properties 
    if (present(printUpdate)) then
      if(printUpdate .eqv. .True.) then
        print *, "  matEnergy at end of timestep =  ", self % matEnergy
        print *, "  T_new =                         ", self % T
      end if
    end if

  end subroutine updateMatIMC

  !!
  !! Material update for ISMC calculation
  !!
  subroutine updateMatISMC(self, tallyEnergy, printUpdate)
    class(baseMgIMCMaterial), intent(inout) :: self
    real(defReal), intent(in)               :: tallyEnergy
    real(defReal)                           :: beta, zeta
    logical(defBool), intent(in), optional  :: printUpdate

    ! Update material internal energy
    self % matEnergy = tallyEnergy

    !if(self % matEnergy <= 0.3) self % matEnergy = 0.3

    ! Update material temperature
    self % T = self % tempFromEnergy()

    ! Update ISMC equivalent of fleck factor
    beta = 4*radiationConstant * self % T**3 / poly_eval(self % cv, self % T)
    self % eta  =   radiationConstant * self % T**4 / self % matEnergy
    zeta = beta - self % eta
    self % fleck = 1 / (1 + zeta*self % sigmaP*lightSpeed*self % deltaT)

    ! Update sigma
    call self % sigmaFromTemp()

    ! Print updated properties 
    if (present(printUpdate)) then
      if(printUpdate .eqv. .True.) then
        print *, "  matEnergy at end of timestep =  ", self % matEnergy
        print *, "  T_new =                         ", self % T
      end if
    end if

  end subroutine updateMatISMC

  !!
  !! Calculate the temperature of material from internal energy
  !!
  function tempFromEnergy(self) result(T)
    class(baseMgIMCMaterial), intent(inout) :: self
    real(defReal)                           :: T, energyDens

    energyDens = self % matEnergy / self % V
    T = poly_solve(self % updateEqn, self % cv, self % T, energyDens)

  end function tempFromEnergy

  !!
  !! Calculate sigma from current temp
  !!
  subroutine sigmaFromTemp(self)
    class(baseMgIMCMaterial), intent(inout) :: self

    self % sigmaP = poly_eval(self % planckEqn, self % T)

    self % data(CAPTURE_XS,:) = poly_eval(self % absEqn, self % T)
    self % data(IESCATTER_XS,:) = poly_eval(self % scattEqn, self % T)
    self % data(TOTAL_XS,:) = self % data(CAPTURE_XS,:) + self % data(IESCATTER_XS,:)
    self % data(PLANCK_XS,:) = poly_eval(self % planckEqn, self % T)

  end subroutine sigmaFromTemp


  !!
  !! Return the energy to be emitted during time step, E_r
  !!
  function getEmittedRad(self) result(emittedRad)
    class(baseMgIMCMaterial), intent(inout) :: self
    real(defReal)                           :: U_r, emittedRad

    U_r = radiationConstant * (self % T)**4

    emittedRad = lightSpeed * self % deltaT * self % sigmaP * self % fleck * U_r * self % V

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
  function getEnergyDens(self) result(energyDens)
    class(baseMgIMCMaterial), intent(inout) :: self
    real(defReal)                           :: energyDens

    energyDens = poly_eval(self % updateEqn, self % T)

  end function getEnergyDens

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
  subroutine setType(self, calcType)
    class(baseMgIMCMaterial), intent(inout) :: self
    integer(shortInt), intent(in)           :: calcType
    real(defReal)                           :: beta, zeta
    character(100), parameter               :: Here = 'setType (baseMgIMCMaterial_class.f90)'

    if(calcType /= IMC .and. calcType /= ISMC) call fatalError(Here, 'Invalid calculation type')

    self % calcType = calcType

    ! If ISMC, recalculate Fleck
    if(self % calcType == ISMC) then
      beta = 4*radiationConstant * self % T**3 / poly_eval(self % cv, self % T)
      self % eta = radiationConstant * self % T**4 / self % matEnergy
      zeta = beta - self % eta
      self % fleck = 1 / (1 + zeta*self % sigmaP*lightSpeed*self % deltaT)
    end if

  end subroutine setType

end module baseMgIMCMaterial_class
