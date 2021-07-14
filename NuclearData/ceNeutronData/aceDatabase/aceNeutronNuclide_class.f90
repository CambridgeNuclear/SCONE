module aceNeutronNuclide_class

  use numPrecision
  use endfConstants
  use universalVariables
  use genericProcedures, only : fatalError, numToChar, binarySearch
  use RNG_class,         only : RNG
  use aceCard_class,     only : aceCard
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
  !!   ZAID           -> ZZAAA.TTc ID of the ACE card of the nuclide
  !!   eGrid          -> Energy grid for the XSs
  !!   mainData       -> Array of XSs that are required in ceNeutronMicroXSs, that is
  !!     (total, capture, escatter, iescatter, fission, nuFission)
  !!   MTdata         -> array of 'reactionMT's with data for all MT reactions in the nuclide
  !!     only reactions 1:nMT are active, that is can be sampled during tracking
  !!   nMT            -> number of active MT reactions that produce 2nd-ary neutrons
  !!   idxMT          -> intMap that maps MT -> index in MTdata array
  !!   elasticScatter -> reactionHandle with data for elastic scattering
  !!   fission        -> reactionHandle with fission data (may be uninitialised)
  !!   urrE           -> energy boundaries of probability tables. It's zero if tables are off
  !!   probTab        -> probability tables for ures
  !!   hasProbTab     -> probability tables flag, it's false by default
  !!   IFF            -> ures probability tables multiplication factor flag
  !!
  !! Interface:
  !!   ceNeutronNuclide Interface
  !!   search    -> search energy grid and return index and interpolation factor
  !!   totalXS   -> return totalXS given index and interpolation factor
  !!   microXSs  -> return interpolated ceNeutronMicroXSs package given index and inter. factor
  !!   getUrrXSs -> return ceNeutronMicroXSs accounting for ures probability tables
  !!   init      -> build nuclide from aceCard
  !!   init_urr  -> build list and mapping of nuclides to maintain temperature correlation
  !!                when reading ures probability tables
  !!   display   -> print information about the nuclide to the console
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
    real(defReal), dimension(2) :: urrE
    type(urrProbabilityTables)  :: probTab
    logical(defBool)            :: hasProbTab = .false.
    integer(shortInt)           :: IFF = 0


  contains
    ! Superclass Interface
    procedure :: invertInelastic
    procedure :: xsOf
    procedure :: kill

    ! Local interface
    procedure :: search
    procedure :: totalXS
    procedure :: microXSs
    procedure :: getUrrXSs
    procedure :: init
    procedure :: init_urr
    procedure :: display

  end type aceNeutronNuclide

contains

  !!
  !! Invert PDF of inelastic stattering
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

    ! Obtain bin index and interpolation factor
    if (nuclideCache(self % getNucIdx()) % E_tot == E) then
      idx = nuclideCache(self % getNucIdx()) % idx
      f   = nuclideCache(self % getNucIdx()) % f

    else
      call self % search(idx, f, E)

    end if

    ! Get inelastic XS
    XS = self % mainData(IESCATTER_XS, idx+1) * f + (ONE-f) * self % mainData(IESCATTER_XS, idx)

    ! Invert
    XS = XS * rand % get()
    do i=1,self % nMT
      ! Get index in MT reaction grid
      idxT = idx - self % MTdata(i) % firstIdx + 1
      if( idxT < 0 ) cycle

      ! Get top and bottom XS
      topXS = self % MTdata(i) % xs(idxT+1)
      bottomXS = self % MTdata(i) % xs(idxT)

      ! Decrement total inelastic and exit if sampling is finished
      XS = XS - topXS * f - (ONE-f) * bottomXS
      if(XS <= ZERO) then
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
  !! TODO: This is quite rought implementation. Improve it!
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

      if(self % isFissile()) then
        xss % fission   = data(FISSION_XS, 2) * f + (ONE-f) * data(FISSION_XS, 1)
        xss % nuFission = data(NU_FISSION, 2) * f + (ONE-f) * data(NU_FISSION, 1)
      else
        xss % fission   = ZERO
        xss % nuFission = ZERO
      end if
    end associate

  end subroutine microXSs

  !!
  !! Return neutronMicroXSs read from probability tables.
  !!
  !! NOTE: The IOA flag for other absorption cross sections is ignored.
  !!       The total xs is not read from tables, but calculated from the other xss.
  !!
  !! Does not perform any check for valid input!
  !!
  !! Args:
  !!   E [in]      -> Energy of ingoing neutron
  !!   xi [in]     -> Random number
  !!   xss [inout] -> XSs package to store interpolated values
  !!
  subroutine getUrrXSs(self, E, xi, xss)
    class(aceNeutronNuclide), intent(in) :: self
    real(defReal), intent(in)            :: E
    real(defReal), intent(in)            :: xi
    type(neutronMicroXSs), intent(inout) :: xss
    real(defReal), dimension(3)          :: val

    ! Read tables
    call self % probTab % sampleXSs(E, xi, val)

    ! Check if multiplication factor is on
    if (self % IFF == 1) then
      val(1) = xss % elasticScatter * val(1)
      val(2) = xss % capture * val(2)
      val(3) = xss % fission * val(3)
    end if

    ! Check inelastic scattering flag
    if (self % probTab % ILF < 0) xss % inelasticScatter = ZERO

    xss % elasticScatter   = val(1)
    xss % capture          = val(2)

    if(self % isFissile()) then
      xss % nuFission = xss % nuFission/xss % fission * val(3)
      xss % fission   = val(3)
    else
      xss % fission   = ZERO
      xss % nuFission = ZERO
    end if

    ! Calculate total from the other values
    xss % total = xss % elasticScatter + xss % inelasticScatter + xss % capture + &
                  xss % fission

  end subroutine getUrrXSs

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

    ! Create a stack of MT reactions, devide them into ones that produce 2nd-ary
    ! particlues and pure absorbtion
    associate (MTs => ACE % getScatterMTs())
      do i=1,size(MTs)
        if (MTs(i) == N_ANYTHING) cycle
        call scatterMT % push(MTs(i))
      end do
    end associate

    associate (MTs => [ACE % getFissionMTs(), ACE % getCaptureMTs()])
      do i=1,size(MTs)
        if(MTs(i) == N_FISSION) cycle ! MT=18 is already included with FIS block
        call absMT % push(MTs(i))
      end do
    end associate

    ! Allocate space
    allocate(self % MTdata(scatterMT % size() + absMT % size()))

    ! Load scattering reactions
    N = scatterMT % size()
    self % nMT = N
    do i =1,N
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
    do i=1,self % nMT
      do j=1,size(self % mainData, 2)
        ! Find bottom and Top of the grid
        bottom = self % MTdata(i) % firstIdx
        top    = size(self % MTdata(i) % xs)
        if( j>= bottom .and. j <= top + bottom) then
          self % mainData(IESCATTER_XS, j) = self % mainData(IESCATTER_XS, j) + &
                                             self % MTdata(i) % xs(j-bottom + 1)
        end if
      end do
    end do

    ! Recalculate totalXS
    if(self % isFissile()) then
      K = FISSION_XS
    else
      K = CAPTURE_XS
    end if
    self % mainData(TOTAL_XS, :) = sum(self % mainData(ESCATTER_XS:K,:),1)

    ! Load Map of MT -> local index of a reaction
    do i=1,size(self % MTdata)
      call self % idxMT % add(self % MTdata(i) % MT, i)
    end do

    ! Include main reaction (in mainData) as -ve entries
    call self % idxMT % add(N_TOTAL,-TOTAL_XS)
    call self % idxMT % add(N_N_ELASTIC, -ESCATTER_XS)
    call self % idxMT % add(N_DISAP, -CAPTURE_XS)

    if(self % isFissile()) then
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
  subroutine init_urr(self, ACE)
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

  end subroutine init_urr

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
