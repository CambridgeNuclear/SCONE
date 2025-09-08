module aceNeutronNuclide_class

  use numPrecision
  use endfConstants
  use universalVariables
  use genericProcedures, only : fatalError, numToChar, binarySearch
  use RNG_class,         only : RNG
  use aceCard_class,     only : aceCard
  use aceSabCard_class,  only : aceSabCard
  use stack_class,       only : stackInt
  use intMap_class,      only : intMap

  ! Nuclear Data Interfaces
  use ceNeutronDatabase_inter,      only : ceNeutronDatabase
  use nuclideHandle_inter,          only : nuclideHandle
  use ceNeutronNuclide_inter,       only : ceNeutronNuclide, kill_super => kill
  use uncorrelatedReactionCE_inter, only : uncorrelatedReactionCE
  use elasticNeutronScatter_class,  only : elasticNeutronScatter
  use fissionCE_class,              only : fissionCE
  use neutronScatter_class,         only : neutronScatter
  use pureCapture_class,            only : pureCapture
  use neutronXSPackages_class,      only : neutronMicroXSs

  ! Unresolved resonances treatment
  use urrProbabilityTables_class,   only : urrProbabilityTables

  ! S(alpha,beta) data
  use thermalScatteringData_class,  only : thermalData

  ! CE NEUTRON CACHE
  use ceNeutronCache_mod,           only : nuclideCache

  implicit none
  private

  !!
  !! Public Pointer Cast
  !!
  public aceNeutronNuclide_CptrCast
  public aceNeutronNuclide_TptrCast


  ! Grid location parameters
  integer(shortInt), parameter :: TOTAL_XS      = 1
  integer(shortInt), parameter :: ESCATTER_XS   = 2
  integer(shortInt), parameter :: IESCATTER_XS  = 3
  integer(shortInt), parameter :: CAPTURE_XS    = 4
  integer(shortInt), parameter :: FISSION_XS    = 5
  integer(shortInt), parameter :: NU_FISSION    = 6


  !!
  !! Groups data related to an MT reaction
  !!
  !! Public Members:
  !!   MT         -> MT number of the reaction
  !!   firstIdx   -> first index occupied by the data on nuclide energy grid
  !!   xs         -> array of reaction XSs
  !!   kinematics -> reaction object for this MT reajction. Contains data about 2nd-ary
  !!     particle distributions
  !!
  type, public :: reactionMT
    integer(shortInt)                         :: MT       = 0
    integer(shortInt)                         :: firstIdx = 0
    real(defReal),dimension(:),allocatable    :: xs
    class(uncorrelatedReactionCE),allocatable :: kinematics
  end type reactionMT

  !!
  !! CE Neutron Nuclide Implementation that follows directly from the specification of ACE data
  !!
  !! NOTE:
  !!   IGNORES MT=5 (N_ANYTHING) in its current implementation !
  !!   IN JEF 3.1.1 It will Reduce accuracy of Tc99 collision processing
  !!
  !! Public Members:
  !!   ZAID             -> ZZAAA.TTc ID of the ACE card of the nuclide
  !!   eGrid            -> Energy grid for the XSs
  !!   mainData         -> Array of XSs that are required in ceNeutronMicroXSs, that is
  !!     (total, capture, escatter, iescatter, fission, nuFission)
  !!   MTdata           -> array of 'reactionMT's with data for all MT reactions in the nuclide
  !!     only reactions 1:nMT are active, that is can be sampled during tracking
  !!   nMT              -> number of active MT reactions that produce 2nd-ary neutrons
  !!   idxMT            -> intMap that maps MT -> index in MTdata array
  !!   elasticScatter   -> reactionHandle with data for elastic scattering
  !!   fission          -> reactionHandle with fission data (may be uninitialised)
  !!   urrE             -> energy boundaries of probability tables. It's zero if tables are off
  !!   probTab          -> probability tables for ures
  !!   hasProbTab       -> probability tables flag, it's false by default
  !!   IFF              -> ures probability tables multiplication factor flag
  !!   hasThData        -> thermal scattering flag, it's false by default
  !!   thData           -> S(a,b) thermal data array to store XSs and outgoing distributions
  !!   stochasticMixing -> flag to indicate whether S(a,b) stochastic interpolation is performed
  !!   SabEl            -> energy boundaries of elastic S(a,b) data
  !!   SabInel          -> energy boundaries of inelastic S(a,b) data
  !!
  !! Interface:
  !!   ceNeutronNuclide Interface
  !!   search           -> search energy grid and return index and interpolation factor
  !!   totalXS          -> return totalXS given index and interpolation factor
  !!   scatterXS        -> return elastic scattering XS given index and interpolation factor
  !!   microXSs         -> return interpolated ceNeutronMicroXSs package given index and inter. factor
  !!   getUrrXSs        -> return ceNeutronMicroXSs accounting for ures probability tables
  !!   getThXSs         -> return ceNeutronMicroXSs accounting for S(a,b) scattering treatment
  !!   getMajXS         -> returns a majorant cross section on request within an energy range given as input
  !!   init             -> build nuclide from aceCard
  !!   initUrr          -> build list and mapping of nuclides to maintain temperature correlation
  !!                       when reading ures probability tables
  !!   initSab          -> builds S(a,b) properties from aceCard
  !!   returnSabPointer -> returns pointer to the appropriate S(a,b) data during stochastic mixing
  !!   display          -> print information about the nuclide to the console
  !!
  type, public, extends(ceNeutronNuclide) :: aceNeutronNuclide
    character(nameLen)                          :: ZAID    = ''
    real(defReal), dimension(:), allocatable    :: eGrid
    real(defReal), dimension(:,:), allocatable  :: mainData
    type(reactionMT), dimension(:), allocatable :: MTdata
    integer(shortInt)                           :: nMT     = 0
    type(intMap)                                :: idxMT

    type(elasticNeutronScatter) :: elasticScatter
    type(fissionCE)             :: fission

    ! URR probability tables
    real(defReal), dimension(2) :: urrE = ZERO
    type(urrProbabilityTables)  :: probTab
    logical(defBool)            :: hasProbTab = .false.
    integer(shortInt)           :: IFF = 0

    ! S(alpha,beta)
    logical(defBool)            :: hasThData = .false.
    logical(defBool)            :: stochasticMixing = .false.
    real(defReal), dimension(2) :: SabEl = ZERO
    real(defReal), dimension(2) :: SabInel = ZERO
    type(thermalData), dimension(:), allocatable :: thData

  contains
    ! Superclass Interface
    procedure :: invertInelastic
    procedure :: xsOf
    procedure :: elScatteringXS
    procedure :: needsSabEl
    procedure :: kill

    ! Local interface
    procedure :: search
    procedure :: totalXS
    procedure :: scatterXS
    procedure :: microXSs
    procedure :: getUrrXSs
    procedure :: getThXSs
    procedure :: getMajXS
    procedure :: needsUrr
    procedure :: needsSabInel
    procedure :: init
    procedure :: initUrr
    procedure :: initSab
    procedure :: getSabPointer
    procedure :: getSabTBounds
    procedure :: display

  end type aceNeutronNuclide

contains

  !!
  !! Invert PDF of inelastic stattering
  !! NOTE: if S(a,b) thermal scattering treatment is on, it is used in place of
  !!       all the other scattering reactions
  !!
  !! TODO: This is quite rought implementation. Improve it!
  !!
  !! See ceNeutronNuclide documentation
  !!
  function invertInelastic(self, E, rand) result(MT)
    class(aceNeutronNuclide), intent(in) :: self
    real(defReal), intent(in)            :: E
    class(RNG), intent(inout)            :: rand
    integer(shortInt)                    :: MT
    integer(shortInt)                    :: idx, i, idxT
    real(defReal)                        :: f, XS, topXS, bottomXS
    character(100), parameter :: Here = 'invertInelastic (aceNeutronNuclide_class.f90)'

    ! Check if it's thermal inelastic scattering or not
    if (self % needsSabInel(E)) then
      MT = N_N_ThermINEL
      return
    end if

    ! Normal (without S(a,b)) inelastic scattering
    ! Obtain bin index and interpolation factor
    call self % search(idx, f, E)

    ! Get inelastic XS
    XS = self % mainData(IESCATTER_XS, idx+1) * f + (ONE-f) * self % mainData(IESCATTER_XS, idx)

    ! Invert
    XS = XS * rand % get()
    do i = 1,self % nMT
      ! Get index in MT reaction grid
      idxT = idx - self % MTdata(i) % firstIdx + 1
      if ( idxT < 1 ) cycle

      ! Get top and bottom XS
      topXS = self % MTdata(i) % xs(idxT+1)
      bottomXS = self % MTdata(i) % xs(idxT)

      ! Decrement total inelastic and exit if sampling is finished
      XS = XS - topXS * f - (ONE-f) * bottomXS
      if (XS <= ZERO) then
        MT = self % MTdata(i) % MT
        return
      end if
    end do

    ! Sampling must have failed. Throw fatalError
    call fatalError(Here,'Failed to invert inelastic scattering in nuclide '//trim(self % ZAID))
    MT = 0

  end function invertInelastic

  !!
  !! Return Cross-Section of reaction MT at energy E
  !!
  !! Needs to use ceNeutronCache
  !!
  !! TODO: This is quite rough implementation. Improve it!
  !!
  !! See ceNeutronNuclide documentation
  !!
  function xsOf(self, MT, E) result(xs)
    class(aceNeutronNuclide), intent(in) :: self
    integer(shortInt), intent(in)        :: MT
    real(defReal), intent(in)            :: E
    real(defReal)                        :: xs
    integer(shortInt)                    :: idx, idxMT
    real(defReal)                        :: f, topXS, bottomXS
    character(100), parameter :: Here = 'xsOf (aceNeutronNuclide_class.f90)'

    ! Find the index of MT reaction in nuclide
    idxMT = self % idxMT % getOrDefault(MT, 0)

    ! Error message if not found
    if (idxMT == 0) then
      call fatalError(Here, 'Requested MT: '//numToChar(MT)// &
                             ' is not present in nuclide '//trim(self % ZAID))
    end if

    ! Obtain bin index and interpolation factor
    if (nuclideCache(self % getNucIdx()) % E_tot == E) then
      idx = nuclideCache(self % getNucIdx()) % idx
      f   = nuclideCache(self % getNucIdx()) % f

    else
      call self % search(idx, f, E)

    end if

    ! Obtain value
    if (idxMT > 0) then
      idx = idx - self % MTdata(idxMT) % firstIdx
      if (idx < 0) then
        topXS = ZERO
        bottomXS = ZERO
      else
        topXS    = self % MTdata(idxMT) % xs(idx+1)
        bottomXS = self % MTdata(idxMT) % xs(idx)
      end if

    else
      idxMT = -idxMT
      topXS    = self % mainData(idxMT, idx+1)
      bottomXS = self % mainData(idxMT, idx)
    end if

    xs = topXS * f + (1-f) * bottomXS

  end function xsOf

  !!
  !! Return value of the elastic scattering XS given neutron energy
  !!
  !! See ceNeutronNuclide documentation
  !!
  function elScatteringXS(self, E) result(xs)
    class(aceNeutronNuclide), intent(in) :: self
    real(defReal), intent(in)            :: E
    real(defReal)                        :: xs
    integer(shortInt)                    :: idx
    real(defReal)                        :: f

    ! Find energy index
    call self % search(idx, f, E)
    ! Retrieve cross section
    xs = self % scatterXS(idx, f)

  end function elScatteringXS

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(aceNeutronNuclide), intent(inout) :: self
    integer(shortInt)                       :: i

    ! Call superclass
    call kill_super(self)

    ! Release reactions
    call self % elasticScatter % kill()
    call self % fission % kill()

    if(allocated(self % MTdata)) then
      do i=1,size(self % MTdata)
        call self % MTdata(i) % kinematics % kill()
      end do
    end if

    ! Local killing
    self % ZAID = ''
    self % nMT  = 0
    if(allocated(self % MTdata))   deallocate(self % MTdata)
    if(allocated(self % mainData)) deallocate(self % mainData)
    if(allocated(self % eGrid))    deallocate(self % eGrid)
    call self % idxMT % kill()

  end subroutine kill

  !!
  !! Search energy for grid and interpolation factor for energy E
  !!
  !! Interpolation factor definition:
  !!   f = (E - E_low) / (E_top - E_low)
  !!   E = E_top * f + E_low * (1-f)
  !!
  !! Args:
  !!   idx [out] -> index of the bottom bin for energy E
  !!   f   [out] -> value of the interpolation factor for energy E
  !!   E   [in]  -> Energy to search for [MeV]
  !!
  !! Errors:
  !!   If energy E is beyond range terminate with fatalError
  !!
  subroutine search(self, idx, f, E)
    class(aceNeutronNuclide), intent(in) :: self
    integer(shortInt), intent(out)       :: idx
    real(defReal), intent(out)           :: f
    real(defReal), intent(in)            :: E
    character(100), parameter :: Here = 'search (aceNeutronNuclide_class.f90)'

    idx = binarySearch(self % eGrid, E)
    if(idx <= 0) then
      call fatalError(Here,'Failed to find energy: '//numToChar(E)//&
                           ' for nuclide '// trim(self % ZAID))
    end if

    associate(E_top => self % eGrid(idx + 1), E_low  => self % eGrid(idx))
      f = (E - E_low) / (E_top - E_low)
    end associate

  end subroutine search

  !!
  !! Return value of the total XS given interpolation factor and index
  !!
  !! Does not perform any check for valid input!
  !!
  !! Args:
  !!   idx [in] -> index of the bottom bin in nuclide Energy-Grid
  !!   f [in]   -> interpolation factor in [0;1]
  !!
  !! Result:
  !!   xs = sigma(idx+1) * f + (1-f) * sigma(idx)
  !!
  !! Errors:
  !!   Invalid idx beyond array bounds -> undefined behaviour
  !!   Invalid f (outside [0;1]) -> incorrect value of XS
  !!
  elemental function totalXS(self, idx, f) result(xs)
    class(aceNeutronNuclide), intent(in) :: self
    integer(shortInt), intent(in)        :: idx
    real(defReal), intent(in)            :: f
    real(defReal)                        :: xs

    xs = self % mainData(TOTAL_XS, idx+1) * f + (ONE-f) * self % mainData(TOTAL_XS, idx)

  end function totalXS

  !!
  !! Return value of the elastic scattering XS given interpolation factor and index
  !!
  !! Does not perform any check for valid input!
  !!
  !! Args:
  !!   idx [in] -> index of the bottom bin in nuclide Energy-Grid
  !!   f [in]   -> interpolation factor in [0;1]
  !!
  !! Result:
  !!   xs = sigma(idx+1) * f + (1-f) * sigma(idx)
  !!
  !! Errors:
  !!   Invalid idx beyond array bounds -> undefined behaviour
  !!   Invalid f (outside [0;1]) -> incorrect value of XS
  !!
  elemental function scatterXS(self, idx, f) result(xs)
    class(aceNeutronNuclide), intent(in) :: self
    integer(shortInt), intent(in)        :: idx
    real(defReal), intent(in)            :: f
    real(defReal)                        :: xs

    xs = self % mainData(ESCATTER_XS, idx+1) * f + (ONE-f) * self % mainData(ESCATTER_XS, idx)

  end function scatterXS

  !!
  !! Return interpolated neutronMicroXSs package for the given interpolation factor and index
  !!
  !! Does not perform any check for valid input!
  !!
  !! Args:
  !!   xss [out] -> XSs package to store interpolated values
  !!   idx [in]  -> index of the bottom bin in nuclide Energy-Grid
  !!   f [in]    -> interpolation factor in [0;1]
  !!
  !! Errors:
  !!   Invalid idx beyond array bounds -> undefined behaviour
  !!   Invalid f (outside [0;1]) -> incorrect value of XSs
  !!
  elemental subroutine microXSs(self, xss, idx, f)
    class(aceNeutronNuclide), intent(in) :: self
    type(neutronMicroXSs), intent(out)   :: xss
    integer(shortInt), intent(in)        :: idx
    real(defReal), intent(in)            :: f

    associate (data => self % mainData(:,idx:idx+1))

      xss % total            = data(TOTAL_XS, 2)  * f + (ONE-f) * data(TOTAL_XS, 1)
      xss % elasticScatter   = data(ESCATTER_XS, 2)  * f + (ONE-f) * data(ESCATTER_XS, 1)
      xss % inelasticScatter = data(IESCATTER_XS, 2) * f + (ONE-f) * data(IESCATTER_XS, 1)
      xss % capture          = data(CAPTURE_XS, 2)   * f + (ONE-f) * data(CAPTURE_XS, 1)

      if (self % isFissile()) then
        xss % fission   = data(FISSION_XS, 2) * f + (ONE-f) * data(FISSION_XS, 1)
        xss % nuFission = data(NU_FISSION, 2) * f + (ONE-f) * data(NU_FISSION, 1)
      else
        xss % fission   = ZERO
        xss % nuFission = ZERO
      end if
    end associate

  end subroutine microXSs

  !!
  !! Return interpolated neutronMicroXSs package for the given interpolation factor and index
  !! including thermal scattering data
  !!
  !! Does not perform any check for valid input!
  !!
  !! NOTE: It recalculates the total cross section given the partials
  !!
  !! Args:
  !!   xss [out]    -> XSs package to store interpolated values
  !!   idx [in]     -> index of the bottom bin in nuclide Energy-Grid
  !!   f [in]       -> interpolation factor in [0;1]
  !!   E [in]       -> Energy of ingoing neutron
  !!   kT [in]      -> Local material thermal energy
  !!   rand [inout] -> RNG for stochastic mixing
  !!
  !! Errors:
  !!   Invalid idx beyond array bounds -> undefined behaviour
  !!   Invalid f (outside [0;1]) -> incorrect value of XSs
  !!
  subroutine getThXSs(self, xss, idx, f, E, kT, rand)
    class(aceNeutronNuclide), intent(in) :: self
    type(neutronMicroXSs), intent(out)   :: xss
    integer(shortInt), intent(in)        :: idx
    real(defReal), intent(in)            :: f
    real(defReal), intent(in)            :: E
    class(RNG), intent(inout)            :: rand
    real(defReal), intent(in)            :: kT
    type(thermalData), pointer           :: sabPtr
    integer(shortInt)                    :: sabIdx

    associate (data => self % mainData(:,idx:idx+1))

      ! Retrieve capture and fission cross sections as usual
      xss % capture = data(CAPTURE_XS, 2) * f + (ONE-f) * data(CAPTURE_XS, 1)

      if (self % isFissile()) then
        xss % fission   = data(FISSION_XS, 2) * f + (ONE-f) * data(FISSION_XS, 1)
        xss % nuFission = data(NU_FISSION, 2) * f + (ONE-f) * data(NU_FISSION, 1)
      else
        xss % fission   = ZERO
        xss % nuFission = ZERO
      end if

      ! Read S(a,b) tables for elastic scatter: return zero if elastic scatter is off.
      ! Default to low temperature without stochastic mixing.
      ! IMPORTANT
      ! The choice of data should be stored somewhere for consistent handling of 
      ! angular distributions, e.g., a cache
      call self % getSabPointer(kT, rand, sabPtr, sabIdx)
      nuclideCache(self % getNucIdx()) % sabIdx = sabIdx

      ! Read S(a,b) tables for elastic scatter: return zero if elastic scatter is off
      xss % elasticScatter = sabPtr % getElXS(E)

      ! If inelastic scatter is on, reads S(a,b) tables for inelastic scatter
      if (self % needsSabInel(E)) then
        xss % inelasticScatter = sabPtr % getInelXS(E)
      else
        xss % inelasticScatter = data(IESCATTER_XS, 2) * f + (ONE-f) * data(IESCATTER_XS, 1)
      end if

    end associate

    ! Calculate total cross sections by summing up the partials
    xss % total = xss % elasticScatter + xss % inelasticScatter + xss % capture + &
                  xss % fission

  end subroutine getThXSs

  !!
  !! Return neutronMicroXSs read from probability tables.
  !!
  !! NOTE: The IOA flag, read from the ACE files, indicates how to determine the cross
  !!       section of 'other absorptions', i.e. all absorptions except for fission and (n,gamma).
  !!       Its contribution is expected to be very small, so here it is ignored as done
  !!       in Serpent and OpenMC as well.
  !!
  !! NOTE: The total xs is not read from tables, but calculated from the partial xss.
  !!
  !! Does not perform any check for valid input!
  !!
  !! Args:
  !!   xss [out] -> XSs package to store interpolated values
  !!   idx [in]  -> index of the bottom bin in nuclide Energy-Grid
  !!   f [in]    -> interpolation factor in [0;1]
  !!   E [in]    -> Energy of ingoing neutron
  !!   xi [in]   -> Random number
  !!
  !! Errors:
  !!   Invalid idx beyond array bounds -> undefined behaviour
  !!   Invalid f (outside [0;1]) -> incorrect value of XSs
  !!
  subroutine getUrrXSs(self, xss, idx, f, E, xi)
    class(aceNeutronNuclide), intent(in) :: self
    type(neutronMicroXSs), intent(out)   :: xss
    integer(shortInt), intent(in)        :: idx
    real(defReal), intent(in)            :: f
    real(defReal), intent(in)            :: E
    real(defReal), intent(in)            :: xi
    real(defReal), dimension(3)          :: val

    ! Read table values for elastic scattering, capture and fission
    call self % probTab % sampleXSs(E, xi, val)

    associate (data => self % mainData(:,idx:idx+1))

      ! Retrieve fission related cross sections as usual
      if (self % isFissile()) then
        xss % fission   = data(FISSION_XS, 2) * f + (ONE-f) * data(FISSION_XS, 1)
        xss % nuFission = data(NU_FISSION, 2) * f + (ONE-f) * data(NU_FISSION, 1)
      else
        xss % fission   = ZERO
        xss % nuFission = ZERO
      end if

      ! Check if flag for multiplication factor (IFF) is true, and apply it to elastic scattering,
      ! capture and fission
      if (self % IFF == 1) then
        xss % elasticScatter   = data(ESCATTER_XS, 2)  * f + (ONE-f) * data(ESCATTER_XS, 1)
        xss % capture          = data(CAPTURE_XS, 2)   * f + (ONE-f) * data(CAPTURE_XS, 1)

        val(1) = xss % elasticScatter * val(1)
        val(2) = xss % capture * val(2)
        val(3) = xss % fission * val(3)
      end if

      ! Check the value of the inelastic scattering flag (ILF). The inelastic scattering cross sections
      ! is treated normally if ILF >= 0
      if (self % probTab % ILF < 0) then
        xss % inelasticScatter = ZERO
      else
        xss % inelasticScatter = data(IESCATTER_XS, 2) * f + (ONE-f) * data(IESCATTER_XS, 1)
      end if

    end associate

    ! Update cross section values
    xss % elasticScatter   = val(1)
    xss % capture          = val(2)

    if(self % isFissile()) then
      xss % nuFission = xss % nuFission/xss % fission * val(3)
      xss % fission   = val(3)
    end if

    ! Calculate total cross section from the partial cross sections
    xss % total = xss % elasticScatter + xss % inelasticScatter + xss % capture + &
                  xss % fission

  end subroutine getUrrXSs

  !! Function to calculate the maximum elastic scattering cross section within
  !! an energy range given by an upper and lower energy bound.
  !!
  !! Args:
  !!   eLower [in]  -> Lower bound of energy range
  !!   eUpper [in]  -> Upper bound of energy range
  !!   MT     [in]  -> MT number of the requested reaction cross section
  !!   maj [out]    -> Maximum scattering cross section within energy range
  !!
  function getMajXS(self, eLower, eUpper, MT) result (maj)
    class(aceNeutronNuclide), intent(in)  :: self
    real(defReal), intent(in)             :: eLower
    real(defReal), intent(in)             :: eUpper
    integer(shortInt), intent(in)         :: MT
    real(defReal)                         :: maj
    integer(shortInt)                     :: reaction, idx
    real(defReal)                         :: f, E, xs
    character(100), parameter :: Here = 'getMajXS (aceNeutronNuclide_class.f90)'

    ! Select desired reaction based on requested MT number
    select case (MT)

      case (N_TOTAL)
        reaction = TOTAL_XS

      case (N_N_ELASTIC)
        reaction = ESCATTER_XS

      case default
        call fatalError(Here, 'Unsupported MT number requested: '//numToChar(MT))

    end select

    ! Search for idx, f, and xs for the lower energy limit
    call self % search(idx, f, eLower)

    ! Conservative: choose the xs at the energy point before the lower energy limit
    maj = self % mainData(reaction, idx)

    majorantLoop: do

      ! Increase index
      idx = idx + 1

      ! Find XS and energy at index
      xs = self % mainData(reaction, idx)
      E  = self % eGrid(idx)

      ! Compare cross sections and possibly update majorant
      maj = max(xs, maj)

      ! Exit loop after getting to the upper energy limit
      if (E >= eUpper) exit majorantLoop

    end do majorantLoop

  end function getMajXS

  !!
  !! Function that checks whether this nuclide at the provided energy should
  !! read unresolved resonance probability tables or not
  !!
  !! Args:
  !!   E [in] -> incident neutron energy
  !!
  !! Returns true or false
  !!
  elemental function needsUrr(self, E) result(doesIt)
    class(aceNeutronNuclide), intent(in)  :: self
    real(defReal), intent(in)             :: E
    logical(defBool)                      :: doesIt

    doesIt = self % hasProbTab .and. E >= self % urrE(1) .and. E <= self % urrE(2)

  end function needsUrr

  !!
  !! Function that checks whether or not this nuclide at the provided energy should
  !! have S(a,b) inelastic scattering data
  !!
  !! Args:
  !!   E [in] -> incident neutron energy
  !!
  !! Returns true or false
  !!
  elemental function needsSabInel(self, E) result(doesIt)
    class(aceNeutronNuclide), intent(in)  :: self
    real(defReal), intent(in)             :: E
    logical(defBool)                      :: doesIt

    doesIt = self % hasThData .and. E >= self % SabInel(1) .and. E <= self % SabInel(2)

  end function needsSabInel

  !!
  !! Function that checks whether or not this nuclide at the provided energy should
  !! have S(a,b) elastic scattering data
  !!
  !! See ceNeutronNuclide documentation
  !!
  elemental function needsSabEl(self, E) result(doesIt)
    class(aceNeutronNuclide), intent(in)  :: self
    real(defReal), intent(in)             :: E
    logical(defBool)                      :: doesIt

    doesIt = self % hasThData .and. E >= self % SabEl(1) .and. E <= self % SabEl(2)

  end function needsSabEl

  !!
  !! Initialise from an ACE Card
  !!
  !! Args:
  !!   ACE [inout]   -> ACE card
  !!   nucIdx [in]   -> nucIdx of the nuclide
  !!   database [in] -> pointer to ceNeutronDatabase to set XSs to cache
  !!
  !! Errors:
  !!   FatalError if ACE card has NU data but no fission MTs
  !!
  subroutine init(self, ACE, nucIdx, database)
    class(aceNeutronNuclide), intent(inout)       :: self
    class(aceCard), intent(inout)                 :: ACE
    integer(shortInt), intent(in)                 :: nucIdx
    class(ceNeutronDatabase), pointer, intent(in) :: database
    integer(shortInt)                             :: Ngrid, N, K, i, j, MT, bottom, top
    type(stackInt)                                :: scatterMT, absMT
    character(100), parameter :: Here = "init (aceNeutronNuclide_class.f90)"

    ! Reset nuclide just in case
    call self % kill()

    ! Load ACE ZAID
    self % ZAID = ACE % ZAID

    ! Read key data into the superclass
    call self % set( fissile  = ACE % isFissile(), &
                     mass     = ACE % AW,          &
                     kT       = ACE % TZ,          &
                     nucIdx   = nucIdx,            &
                     database = database )

    ! Get size of the grid
    Ngrid = ACE % gridSize()

    ! Allocate space for main XSs
    if(self % isFissile()) then
      N = 6
    else
      N = 4
    end if
    allocate(self % mainData(N, Ngrid))

    self % mainData = ZERO

    ! Load Main XSs
    self % eGrid =  ACE % ESZ_XS('energyGrid')
    self % mainData(TOTAL_XS,:)     = ACE % ESZ_XS('totalXS')
    self % mainData(ESCATTER_XS,:)  = ACE % ESZ_XS('elasticXS')
    self % mainData(CAPTURE_XS,:)   = ACE % ESZ_XS('absorptionXS')

    ! Get elastic kinematics
    call self % elasticScatter % init(ACE, N_N_ELASTIC)

    ! Load Fission XS data
    ! Set 'bottom' variable to the start index of fission data
    if (self % isFissile()) then

      if (ACE % hasFIS()) then
        ! Generic fission reaction MT=18 is provided
        ! Read XS from the FIS block in ACE card
        N = ACE % firstIdxFiss()
        bottom = N
        K = ACE % numXSPointsFiss()
        self % mainData(FISSION_XS, N:N+K-1) = ACE % xsFiss()

      else
        ! FIS block is missing, so MT=18 is unavaliable
        ! Build fission XS as a sum of MT=19,20,21,38
        ! if they are present in the ACE card
        associate (fissMTs => ACE % getFissionMTs())

          ! Perform sanity check
          if (size(fissMTs) < 1) then
            call fatalError(Here, ACE % ZAID//' seems to have NU data but no fission reactions &
                                  &among its MT numbers')
          end if

          ! Get fission XS via summation
          bottom = Ngrid + 1
          do i = 1, size(fissMTs)
            N = ACE % firstIdxMT(fissMTs(i))
            bottom = min(bottom, N)
            K = ACE % numXsPointsMT(fissMts(i))
            self % mainData(FISSION_XS, N:N+K-1) = self % mainData(FISSION_XS, N:N+K-1) &
                                                 + ACE % xsMT(fissMTs(i))
          end do

        end associate
      end if

      ! Build Fission reaction object
      ! Again select between MT=18 and 19 based on presence of FIS block
      if (ACE % hasFIS()) then
        call self % fission % init(ACE, N_FISSION)
      else
        call self % fission % init(ACE, N_f)
      end if

      ! Calculate nuFission
      do i = bottom, Ngrid
        self % mainData(NU_FISSION,i) = self % mainData(FISSION_XS,i) * &
                                        self % fission % release(self % eGrid(i))
      end do

    end if

    ! Read data for MT reaction

    ! Create a stack of MT reactions, divide them into ones that produce 2nd-ary
    ! particles and pure absorption
    associate (MTs => ACE % getScatterMTs())
      do i = 1,size(MTs)
        if (MTs(i) == N_ANYTHING) cycle
        call scatterMT % push(MTs(i))
      end do
    end associate

    associate (MTs => [ACE % getFissionMTs(), ACE % getCaptureMTs()])
      do i = 1,size(MTs)
        if(MTs(i) == N_FISSION) cycle ! MT=18 is already included with FIS block
        call absMT % push(MTs(i))
      end do
    end associate

    ! Allocate space
    allocate(self % MTdata(scatterMT % size() + absMT % size()))

    ! Load scattering reactions
    N = scatterMT % size()
    self % nMT = N
    do i = 1,N
      call scatterMT % pop(MT)
      self % MTdata(i) % MT       = MT
      self % MTdata(i) % firstIdx = ACE % firstIdxMT(MT)
      self % MTdata(i) % xs       = ACE % xsMT(MT)

      allocate(neutronScatter :: self % MTdata(i) % kinematics)
      call self % MTdata(i) % kinematics % init(ACE, MT)
    end do

    ! Load capture reactions
    K = absMT % size()
    do i = N+1,N+K
      call absMT % pop(MT)
      self % MTdata(i) % MT       = MT
      self % MTdata(i) % firstIdx = ACE % firstIdxMT(MT)
      self % MTdata(i) % xs       = ACE % xsMT(MT)

      allocate(pureCapture :: self % MTdata(i) % kinematics)
      call self % MTdata(i) % kinematics % init(ACE, MT)
    end do

    ! Calculate Inelastic scattering XS
    do i = 1,self % nMT
      do j = 1,size(self % mainData, 2)
        ! Find bottom and Top of the grid
        bottom = self % MTdata(i) % firstIdx
        top    = size(self % MTdata(i) % xs)
        if (j>= bottom .and. j <= top + bottom) then
          self % mainData(IESCATTER_XS, j) = self % mainData(IESCATTER_XS, j) + &
                                             self % MTdata(i) % xs(j-bottom + 1)
        end if
      end do
    end do

    ! Recalculate totalXS
    if (self % isFissile()) then
      K = FISSION_XS
    else
      K = CAPTURE_XS
    end if
    self % mainData(TOTAL_XS, :) = sum(self % mainData(ESCATTER_XS:K,:),1)

    ! Load Map of MT -> local index of a reaction
    do i = 1,size(self % MTdata)
      call self % idxMT % add(self % MTdata(i) % MT, i)
    end do

    ! Include main reaction (in mainData) as -ve entries
    call self % idxMT % add(N_TOTAL,-TOTAL_XS)
    call self % idxMT % add(N_N_ELASTIC, -ESCATTER_XS)
    call self % idxMT % add(N_DISAP, -CAPTURE_XS)

    if (self % isFissile()) then
      call self % idxMT % add(N_FISSION, -FISSION_XS)
    end if

    ! TODO: Uncomment after shrinking is implemented in intMap
    !call self % idxMT % shrink()

  end subroutine init

  !!
  !! Initialise probability tables from ACE card
  !!
  !! Switches off tables for this nuclide if tables have some inconsistency
  !!
  !! Args:
  !!   ACE [inout]   -> ACE card
  !!
  subroutine initUrr(self, ACE)
    class(aceNeutronNuclide), intent(inout) :: self
    class(aceCard), intent(inout)           :: ACE

    ! Read in ACE if the nuclide has URR probability tables
    self % hasProbTab = ACE % hasProbTab()

    if (self % hasProbTab) then
      ! Initialise probability tables
      call self % probTab % init(ACE)
      ! Check if probability tables were read correctly

      if (allocated(self % probTab % eGrid)) then
        self % urrE = self % probTab % getEbounds()
        self % IFF = self % probTab % getIFF()

      else
        ! Something went wrong!
        self % hasProbTab = .false.
        self % urrE = ZERO

      end if

    else
      self % urrE = ZERO

    end if

  end subroutine initUrr

  !!
  !! Initialise thermal scattering tables from ACE card
  !!
  !! Args:
  !!   ACE1 [inout]   -> ACE S(a,b) card
  !!   ACE2 [inout]   -> Optional second ACE S(a,b) card
  !!
  !! Errors:
  !!   fatalError if the inelastic scattering S(a,b) energy grid starts at a
  !!   lower energy than the nuclide energy grid
  !!
  subroutine initSab(self, ACE1, ACE2)
    class(aceNeutronNuclide), intent(inout)    :: self
    class(aceSabCard), intent(inout)           :: ACE1
    class(aceSabCard), intent(inout), optional :: ACE2
    real(defReal), dimension(2)                :: EBounds
    real(defReal)                              :: T1, T2
    type(thermalData)                          :: temp
    character(100), parameter :: Here = "initSab (aceNeutronNuclide_class.f90)"

    if (present(ACE2)) then
      allocate(self % thData(2))
    else
      allocate(self % thData(1))
    end if

    ! Initialise S(a,b) class from ACE file
    call self % thData(1) % init(ACE1)
    self % hasThData = .true.

    ! Initialise energy boundaries
    self % SabInel = self % thData(1) % getEBounds('inelastic')
    self % SabEl = self % thData(1) % getEBounds('elastic')
    
    ! Add second S(a,b) file for stochastic mixing
    if (present(ACE2)) then
      
      self % stochasticMixing = .true.
      call self % thData(2) % init(ACE2)
      
      ! Ensure energy bounds are conservative
      EBounds = self % thData(2) % getEBounds('inelastic')
      if (EBounds(1) > self % SabInel(1)) self % SabInel(1) = EBounds(1)
      if (EBounds(2) < self % SabInel(2)) self % SabInel(2) = EBounds(2)
      
      EBounds = self % thData(2) % getEbounds('elastic')
      if (EBounds(1) > self % SabEl(1)) self % SabEl(1) = EBounds(1)
      if (EBounds(2) < self % SabEl(2)) self % SabEl(2) = EBounds(2)

      ! Identify which data is higher temperature and which is lower
      ! 1 should be lower than 2 - swap if necessary
      T1 = self % thData(1) % getTemperature()
      T2 = self % thData(2) % getTemperature()

      if (T1 > T2) then
        temp = self % thData(1)
        self % thData(1) = self % thData(2)
        self % thData(2) = temp
      end if

    end if

    ! Check consistency of energy grid
    if (self % SabInel(1) < self % eGrid(1)) then
      call fatalError(Here, 'S(a,b) low energy boundary is lower than nuclide first energy point')
    end if

  end subroutine initSab

  !!
  !! Return pointer to Sab reaction data.
  !! If stochastic mixing is active, samples which
  !! set of reaction data to point towards.
  !!
  !! Returns the Sab index for later consistent handling
  !! of angular distributions by storing in a cache
  !!
  subroutine getSabPointer(self, kT, rand, ptr, idx)
    class(aceNeutronNuclide), intent(in), target :: self
    real(defReal), intent(in)                    :: kT
    class(RNG), intent(inout)                    :: rand
    type(thermalData), pointer, intent(out)      :: ptr
    integer(shortInt), intent(out)               :: idx
    real(defReal)                                :: kT1, kT2
    character(100), parameter :: Here = "getSabPointer (aceNeutronNuclide_class.f90)"

    if (self % stochasticMixing) then
      kT1 = self % thData(1) % getTemperature()
      kT2 = self % thData(2) % getTemperature()

      if ((kT < kT1) .or. (kT > kT2)) call fatalError(Here,&
              'Requested temperature '//numToChar(kT)//' not in temperature bounds: '//&
              numToChar(kT1)//' and '//numToChar(kT2))

      if ((kT2 - kT)/(kT2 - kT1) > rand % get()) then
        ptr => self % thData(1)
        idx = 1
      else
        ptr => self % thData(2)
        idx = 2
      end if
    else
      ptr => self % thData(1)
      idx = 1
    end if

  end subroutine getSabPointer

  !!
  !! Return the temperature bounds of S(alpha,beta) data
  !! If only one library, both bounds are the same temperature
  !!
  function getSabTBounds(self) result(kT)
    class(aceNeutronNuclide), intent(in)  :: self
    real(defReal), dimension(2)           :: kT

    kT(1) = self % thData(1) % getTemperature()
    if (self % stochasticMixing) then
      kT(2) = self % thData(2) % getTemperature()
    else
      kT(2) = kT(1)
    end if

  end function getSabTBounds

  !!
  !! A Procedure that displays information about the nuclide to the screen
  !!
  !! Prints:
  !!   nucIdx; Size of Energy grid; Maximum and Minumum energy; Mass; Temperature; ACE ZAID;
  !!   MT numbers used in tranposrt (active MTs) and not used in transport (inactiveMTs)
  !!
  !! NOTE: For Now Formatting is horrible. Can be used only for debug
  !!       MT = 18 may appear as inactive MT but is included from FIS block !
  !!
  !! Args:
  !!   None
  !!
  !! Errors:
  !!   None
  !!
  subroutine display(self)
    class(aceNeutronNuclide), intent(in)  :: self
    integer(shortInt)                     :: N, sMT, allMT
    real(defReal)                         :: E_min, E_max, M, kT

    ! Print Most relevant information
    print *, repeat('#',40)
    print '(A)', "Nuclide: " // trim(self % ZAID)

    N = size(self % eGrid)
    E_min = self % eGrid(1)
    E_max = self % eGrid(N)
    M = self % getMass()
    kT = self % getkT()
    print '(A)', "Mass: " //numToChar(M) //" Temperature: "//numToChar(kT) //" [MeV]"
    print '(A)', "Energy grid: " //numToChar(N)// " points from "// &
                  numToChar(E_min)//" to "// numToChar(E_max) // " [MeV] "

    ! Print MT information
    sMT = self % nMT
    allMT = size(self % MTdata)
    print '(A)', "Active MTs: "  // numToChar(self % MTdata(1:sMT) % MT)
    print '(A)', "Inactive MTs: "// numToChar(self % MTdata(sMT+1:allMT) % MT)
    print '(A)', "This implementation ignores MT=5 (N,anything) !"

  end subroutine display

  !!
  !! Cast nuclideHandle pointer to aceNeutronNuclide type pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class nuclideHandle
  !!
  !! Result:
  !!   Null is source is not of aceNeutronNuclide type
  !!   Pointer to source if source is aceNuclearDatabase type
  !!
  pure function aceNeutronNuclide_TptrCast(source) result(ptr)
    class(nuclideHandle), pointer, intent(in) :: source
    type(aceNeutronNuclide), pointer          :: ptr

    select type(source)
      type is(aceNeutronNuclide)
        ptr => source

      class default
        ptr => null()
    end select

  end function aceNeutronNuclide_TptrCast

  !!
  !! Cast nuclideHandle pointer to aceNeutronNuclide class pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class nuclideHandle
  !!
  !! Result:
  !!   Null is source is not of aceNeutronNuclide class
  !!   Pointer to source if source is aceNuclearDatabase class
  !!
  pure function aceNeutronNuclide_CptrCast(source) result(ptr)
    class(nuclideHandle), pointer, intent(in) :: source
    class(aceNeutronNuclide), pointer         :: ptr

    select type(source)
      class is(aceNeutronNuclide)
        ptr => source

      class default
        ptr => null()
    end select

  end function aceNeutronNuclide_CptrCast


end module aceNeutronNuclide_class
