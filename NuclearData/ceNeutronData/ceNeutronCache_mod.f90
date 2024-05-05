!!
!! Global Cache for CE Neutron Nuclear Data
!!
!! Any active nuclear database writes XSs to this module.
!! Any nuclide or material handles read from this module
!!
!! Public Members:
!!   nuclideCache  -> Array of cached data for nuclides
!!   materialCache -> Array of cached data for materials
!!   majorantCache -> Array of cached data for majorant XS
!!   trackingCache -> Array of cached data for tracking XS
!!   zaidCache     -> Array of cached data for ZAID numbers (used for URR probability tables)
!!
!! NOTE:
!!   Cache arrays are deliberatly not targets. This is becouse there should be no pointers to the
!!   cache. Any update call can change energy of any value so it would not be possible that the
!!   energy of XSs pointed by pointers would not change silently.
!!
module ceNeutronCache_mod

  use numPrecision
  use genericProcedures,       only : fatalError, numToChar
  use neutronXsPackages_class, only : neutronMicroXSs, neutronMacroXSs

  implicit none
  private

  !!
  !! Structure that contains cached data for each CE Neutron Material
  !!
  !! Public Members:
  !!   E_tot  -> Energy of the total XS in xss
  !!   E_tail -> Energy of all XSs in xss except total
  !!   f      -> Interpolation factor for the nuclide at energy E_tot
  !!   idx    -> Index on a nuclide grid for energy E_tot
  !!   xss    -> Cached Cross-Section values
  !!
  type, public :: cacheMatDat
    real(defReal)         :: E_tot  = ZERO
    real(defReal)         :: E_tail = ZERO
    real(defReal)         :: f      = ZERO
    integer(shortInt)     :: idx    = 0
    type(neutronMacroXSs) :: xss
  end type cacheMatDat

  !!
  !! Structure that contains cached data for each CE Neutron Nuclide
  !!
  !! Public Members:
  !!   E_tot  -> Energy of the total XS in xss
  !!   E_tail -> Energy of all XSs in xss except total
  !!   f      -> Interpolation factor for the nuclide at energy E_tot
  !!   idx    -> Index on a nuclide grid for energy E_tot
  !!   xss    -> Cached Cross-Sections values
  !!   needsSabInel -> Flag that tells if the nuclide is using thermal inelastic
  !!                   scattering data
  !!   needsSabEl   -> Flag that tells if the nuclide is using thermal elastic
  !!                   scattering data
  !!   needsUrr     -> Flag that tells if the nuclide is using unresolved resonance
  !!                   probability tables
  !!
  type, public :: cacheNucDat
    real(defReal)         :: E_tot  = ZERO
    real(defReal)         :: E_tail = ZERO
    real(defReal)         :: f      = ZERO
    integer(shortInt)     :: idx    = 0
    type(neutronMicroXSs) :: xss
    logical(defBool)      :: needsSabInel = .false.
    logical(defBool)      :: needsSabEl = .false.
    logical(defBool)      :: needsUrr = .false.
  end type cacheNucDat

  !!
  !! Structure that contains a cross section value
  !!
  !! Public Members:
  !!   E  -> energy of the cross section
  !!   xs -> value of the cross section
  !!
  type, public :: cacheSingleXS
    real(defReal) :: E
    real(defReal) :: xs
  end type cacheSingleXS

  !!
  !! Structure that contains a ZAID random number for probability tables
  !!
  !! Public Members:
  !!   E  -> energy of the majorant
  !!   xi -> value of the random number
  !!
  type, public :: cacheZAID
    real(defReal) :: E  = ZERO
    real(defReal) :: xi = ZERO
  end type cacheZAID

  ! MEMBERS OF THE MODULE ARE GIVEN HERE
  type(cacheMatDat), dimension(:), allocatable, public   :: materialCache
  type(cacheNucDat), dimension(:), allocatable, public   :: nuclideCache
  type(cacheSingleXS), dimension(:), allocatable, public :: majorantCache
  type(cacheSingleXS), dimension(:), allocatable, public :: trackingCache
  type(cacheZAID), dimension(:), allocatable, public     :: zaidCache
  !$omp threadprivate(materialCache, nuclideCache, majorantCache, trackingCache, zaidCache)

  ! Public procedures
  public :: init
  public :: kill


contains


  !!
  !! Initialise Cache to given space
  !!
  !! Args:
  !!   nMat [in]  -> Number of materials
  !!   nNuc [in]  -> Number of nuclides
  !!   nMaj [in]  -> Optional. Number of majorant Xss (Default = 1)
  !!   nZaid [in] -> Optional. Number of nuclides with same ZAID
  !!
  !! Errors:
  !!   fatalError if nMat, nNuc or Nmaj is not a +ve value
  !!
  subroutine init(nMat, nNuc, nMaj, nZaid)
    integer(shortInt), intent(in)           :: nMat
    integer(shortInt), intent(in)           :: nNuc
    integer(shortInt), optional, intent(in) :: nMaj
    integer(shortInt), optional, intent(in) :: nZaid
    integer(shortInt)                       :: nLoc
    character(100),parameter :: Here = 'init (ceNeutronCache_mod.f90)'

    ! Make sure memory is clean
    call kill()

    ! Read default value of majorant XSs
    if (present(nMaj)) then
      nLoc = nMaj
    else
      nLoc = 1
    end if

    ! Chack the provided data
    if (nMat < 1) call fatalError(Here,'Number of materials must be +ve! Not: '//numToChar(nMat))
    if (nNuc < 1) call fatalError(Here,'Number of nuclides must be +ve! Not: '//numToChar(nMat))
    if (nLoc < 1) call fatalError(Here,'Number of majorant XSs must be +ve! Not: '//numToChar(nMat))

    ! Allocate space
    ! Need to do in parallel region to allocate each copy
    !$omp parallel
    allocate(materialCache(nMat))
    allocate(nuclideCache(nNuc))
    allocate(majorantCache(nLoc))
    allocate(trackingCache(1))

    if (present(nZaid)) then
      if (nZaid > 0) then
        allocate(zaidCache(nZaid))
      else
        call fatalError(Here,'Number of zaids must be +ve! Not: '//numToChar(nZaid))
      end if
    end if
    !$omp end parallel

  end subroutine init

  !!
  !! Return Cache Module (Singleton) to uninitialised state
  !!
  subroutine kill()
    ! Need to deallocate on all threads
    !$omp parallel
    if (allocated(materialCache)) deallocate (materialCache)
    if (allocated(nuclideCache))  deallocate (nuclideCache)
    if (allocated(majorantCache)) deallocate (majorantCache)
    if (allocated(trackingCache)) deallocate (trackingCache)
    if (allocated(zaidCache))     deallocate (zaidCache)
    !$omp end parallel
  end subroutine kill


end module ceNeutronCache_mod
