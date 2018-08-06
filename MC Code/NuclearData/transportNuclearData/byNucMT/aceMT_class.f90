module aceMT_class

  use numPrecision
  use endfConstants
  use genericProcedures ,         only : openToRead, binarySearch, interpolate, searchError, fatalError
  use aceCard_class,              only : aceCard
  use RNG_class,                  only : RNG

  use indexStorageMT_class,       only : indexStorageMT

  ! Emission data objects
  use emissionENDF_class,         only : emissionENDF
  use releaseLawENDF_inter,       only : releaseLawENDF
  use releaseLawENDFfactory_func, only : new_totalNu

  ! XS data objects
  use xsEnergyPointNoMT_class,    only : xsEnergyPointNoMT
  use materialDataNoMT_class,     only : materialDataNoMT
  use reactionMT_class,           only : reactionMT

  implicit none
  private

  !!
  !! Class that stores nuclide reaction cross-sections read from ACE data card.
  !!
  type, public :: aceMT
    real(defReal),public      :: atomWeight    !! Atomic weight ratio A/m_neutron
    real(defReal),public      :: temp          !! Temperature of nucleide [MeV]
    character(nameLen),public :: zzId          !! Nuclide ZZid ie. ZZZAAA.nnC
    character(10)             :: date          !! Date the data were processed
    character(70)             :: comment       !! Quick! Guess what is it!
    character(10)             :: MATid         !! MAT indentifier (see p. F-9 Appendic F MCNP 4 manual)

    logical(defBool)   :: isFissile = .false.

    type(emissionENDF)                           :: eScatterKinematics
    type(emissionENDF)                           :: fissionKinematics
    class(releaseLawENDF),allocatable            :: nuData

    real(defReal),dimension(:),allocatable       :: energyGrid         !! Energy grid for xs data
    type(xsEnergyPointNoMT),dimension(:),pointer :: xsData => null()
    ! **** ABOVE is EFFECTIVLY SUPERCLASS
    type(indexStorageMT)                      :: lookupTableMT
    type(reactionMT)                          :: elasticScatter
    type(reactionMT),dimension(:),allocatable :: MTreactions
    integer(shortInt)                         :: nMTscatter     !! Number of scattering MT reactions

    ! -> Translate MT to an index in nuclide
    ! -> Data for elastic scattering
    ! -> Data for react of scattering reactions -> number of indexes
    ! -> Data for capture reactions

  contains
    procedure :: init
    procedure :: energyIdxFor
    procedure :: sampleMuEout
    procedure :: releaseAt
    procedure :: isInCMframe

    ! ** NEW procedures
    procedure :: invertScatter ! ***

    procedure,private :: readXS

  end type aceMT

contains

  !!
  !! Inverts scattering reaction given random number, anyScatter XS,
  !! energy grid idx and interpolation factor
  !!
  function invertScatter(self,r,XSany,idx,f) result(MT)
    class(aceMT), intent(in)      :: self
    real(defReal), intent(in)     :: r
    real(defReal), intent(in)     :: XSany
    integer(shortInt), intent(in) :: idx    !! Index on energy grid
    real(defReal), intent(in)     :: f      !! Interpolation factor for the current energy
    integer(shortInt)             :: MT
    real(defReal)                 :: XS
    integer(shortInt)             :: i
    character(100),parameter  :: Here = 'invertScattering (aceMT_class.f90)'

    ! Calculate scaled any scattering XS
    XS = r * XSany

    XS = XS - self % elasticScatter % getXS(idx,f)

    if( XS < ZERO) then
      MT = N_N_elastic
      return
    end if

    do i=1,self % nMTscatter
      XS = XS - self % MTreactions(i) % getXS(idx,f)
      if( XS < ZERO) then
        MT = self % MTreactions(i) % MT
        return
      end if

    end do

    call fatalError(Here,'Scattering failed to invert')

  end function invertScatter

  !!
  !! Searches for the given energy and returns index of the bottom of the interval
  !!
  function energyIdxFor(self,E) result(idx)
    class(aceMT), intent(in)    :: self
    real(defReal), intent(in)   :: E
    integer(shortInt)           :: idx
    character(100), parameter   :: Here = 'searchFor (aceMT_class.f90)'

    idx = binarySearch(self % energyGrid,E)
    call searchError(idx,Here)

  end function energyIdxFor


  !!
  !! Sample deflection angle and emission energy for a given MT number
  !!
  subroutine sampleMuEout(self,mu,E_out,E_in,rand,MT)
    class(aceMT), intent(in)      :: self
    real(defReal), intent(out)    :: mu
    real(defReal), intent(out)    :: E_out
    real(defReal), intent(in)     :: E_in
    class(RNG), intent(inout)     :: rand
    integer(shortInt), intent(in) :: MT
    integer(shortInt)             :: idx
    character(100), parameter     :: Here = 'sampleMuEout (aceMT_class.f90)'

    select case(MT)
      case (anyScatter, N_N_elastic)
        call self % eScatterKinematics % sampleAngleEnergy(mu,E_out,E_in,rand)

      case (anyCapture, N_disap:) ! Any capture and all MT >= 101
        mu    = ONE
        E_out = E_in

      case (anyFission, N_fission)
        call self % fissionKinematics  % sampleAngleEnergy(mu,E_out,E_in,rand)

      case default
        idx = self % lookupTableMT % getIdx(MT)
        call self % MTreactions(idx)  % sampleAngleEnergy(mu,E_out,E_in,rand)
        !call fatalError(Here,'Unknown MT number')

    end select

  end subroutine sampleMuEout

  !!
  !! Function that returns .true. if emission data for a given MT is stored in Center-Of-Mass frame
  !!
  function isInCMframe(self,MT) result (isIt)
    class(aceMT), intent(in)     :: self
    integer(shortInt), intent(in) :: MT
    logical(defBool)              :: isIt
    integer(shortInt)             :: idx
    character(100), parameter     :: Here = 'isInCMframe (aceMT_class.f90)'

    select case(MT)
      case (anyScatter, N_N_elastic)
        isIt = self % eScatterKinematics % isInCMframe()

      case (anyCapture, N_disap:) ! Any capture and all MT >= 101
        isIt = .false.

      case (anyFission, N_fission)
        isIt = self % fissionKinematics % isInCMframe()

      case default
        idx = self % lookupTableMT % getIdx(MT)
        isIt = self % MTreactions(idx) % isInCMframe()

    end select


  end function isInCMframe


  !!
  !! Sample give value of average neutron emission for given energy and MT number
  !!
  function releaseAt(self,E_in,MT) result (nu)
    class(aceMT), intent(in)      :: self
    real(defReal), intent(in)     :: E_in
    integer(shortInt), intent(in) :: MT
    real(defReal)                 :: nu
    character(100), parameter     :: Here = 'releaseAt (aceMT_class.f90)'

    select case(MT)
      case (anyScatter, N_N_elastic)
        nu = ONE

      case (anyCapture, N_disap:) ! Any capture and all MT >= 101
        nu = ZERO

      case (anyFission, N_fission)
        nu = self % nuData % releaseAt(E_in)

      case default
        call fatalError(Here,'Unknown MT number')

    end select

  end function releaseAt

  !!
  !! Load reaction cross-section data from ACE file.
  !!
  subroutine init(self,filePath,line)
    class(aceMT), intent(inout)             :: self
    character(*), intent(in)                :: filePath  !! Path to file with ACE data card
    integer(shortInt), intent(in)           :: line      !! Line at which the ACE data card begins in the file
    type(aceCard)                           :: ACE

    ! Load ACE card
    call ACE % readFromFile(filePath,line)

    ! Copy nuclide header
    self % ZZid       = ACE % ZAID
    self % atomWeight = ACE % AW
    self % temp       = ACE % TZ
    self % date       = ACE % HD
    self % comment    = ACE % HK
    self % MATid      = ACE % HM

    call self % readXS(ACE)

  end subroutine init

  !!
  !! Read cross-sections into type variables from NXS,JXS and XSS tabels
  !! This implementation ignores reaction other then total,elastic scattering, total capture and
  !! fission.
  !!
  subroutine readXS(self,ACE)
    class(aceMT), intent(inout)                  :: self
    type(aceCard), intent(inout)                 :: ACE
    real(defReal), dimension(:), allocatable     :: xsAnyScatter
    real(defReal), dimension(:), allocatable     :: xsCapture
    real(defReal), dimension(:), allocatable     :: xsFission
    integer(shortInt)                            :: i,j, Ngrid
    integer(shortInt), dimension(:), allocatable :: MTs

    if (allocated(self % energyGrid)) deallocate(self % energyGrid)
    if (associated(self % xsData)) deallocate(self % xsData)

    ! Check if nuclide is fissile
    self % isFissile = ACE % isFissile()

    ! Load Energy Grid
    call ACE % ESZblock(self % energyGrid,'energyGrid')

    ! Load cross-sections from ESZ block of ACE data card.
    call ACE % ESZblock(xsCapture,'absorptionXS')  ! Total capture XS

    ! Find size of energy grid
    Ngrid = size( self % energyGrid)

    ! Load fission data if present
    if (self % isFissile) then
      ! First index on energy grid for which fission data is present
      i = ACE % firstIdxFiss()

      ! Number of consecutive eneteries for data
      j = ACE % numXsPointsFiss()

      ! Allocate adn load fission XSs
      allocate(xsFission(Ngrid))
      xsFission(i:i+j-1) = ACE % xsFiss()

    end if

    ! Load data for all MT reactions
    call self % elasticScatter % init(ACE,N_N_elastic)

    MTs = [ACE % getScatterMTs(), ACE % getFissionMTs(), ACE % getCaptureMTs()]

    ! Allocate space for other reactions
    allocate(self % MTreactions(size(MTs)))

    do i=1,size(MTs)
      call self % MTreactions(i) % init(ACE,MTs(i))

    end do

    ! Store MT indeces in lookup table
    call self % lookupTableMT % init(MTs)

    self % nMTscatter = size(ACE % getScatterMTs())

    ! Calculate anyScatter xs
    allocate(xsAnyScatter(size(self % energyGrid)))

    xsAnyScatter = ZERO

    do i =1,self % nMTscatter
      call self % MTreactions(i) % addOnGrid(xsAnyScatter)
    end do

    call self % elasticScatter % addOnGrid(xsAnyScatter)

    ! Allocate xsData
    allocate(self % xsData(size(self % energyGrid)))

    ! Move cross section data into xsEnergyPointType
    call self % xsData % load (xsAnyScatter,xsCapture,xsFission)

    ! Read elastic scattering angles
    call self % eScatterKinematics % init(ACE,N_N_elastic)

    if (self % isFissile) then
      call self % fissionKinematics % init(ACE,N_fission)

      ! Read and allocate NU data
      allocate(self % nuData, source = new_totalNu(ACE))

    end if

  end subroutine readXS


end module aceMT_class



