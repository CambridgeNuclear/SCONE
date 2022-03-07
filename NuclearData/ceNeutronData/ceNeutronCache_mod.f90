!!
!! Global Cache for CE Neutron Nuclear Data
!!
!! Any active nuclear database writes XSs to this module.
!! Any nuclide or material handles read from this module
!!
!! Public Members:
!!   nuclideCache  -> Array of cached data for nuclides
!!   matCache      -> Array of cached data for materials
!!   majorantCache -> Array of cached data for majorant XS
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
  !!   E_uni  -> Energy idx and f refer to, when material-unionised grids are used
  !!   f      -> Interpolation factor for the nuclide at energy E_tot
  !!   idx    -> Index on a nuclide grid for energy E_tot
  !!   xss    -> Cached Cross-Section values
  !!
  type, public :: cacheMatDat
    real(defReal)         :: E_tot  = ZERO
    real(defReal)         :: E_tail = ZERO
    real(defReal)         :: E_uni  = ZERO
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
  !!   E_mat  -> Energy idx and f refer to, when unionised double-indexing grids are used
  !!   f      -> Interpolation factor for the nuclide at energy E_tot
  !!   idx    -> Index on a nuclide grid for energy E_tot
  !!   xss    -> Cached Cross-Sections values
  !!
  type, public :: cacheNucDat
    real(defReal)         :: E_tot  = ZERO
    real(defReal)         :: E_tail = ZERO
    real(defReal)         :: E_mat  = ZERO
    real(defReal)         :: f      = ZERO
    integer(shortInt)     :: idx    = 0
    type(neutronMicroXSs) :: xss
  end type cacheNucDat

  !!
  !! Structure that contains a Majorant XS
  !!
  !! Public Members:
  !!   E  -> energy of the majorant
  !!   xs -> value of the majorant
  !!
  type, public :: cacheMajorant
    real(defReal) :: E
    real(defReal) :: xs
  end type cacheMajorant

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
  type(cacheMajorant), dimension(:), allocatable, public :: majorantCache
  type(cacheZAID), dimension(:), allocatable, public     :: zaidCache
  !$omp threadprivate(materialCache, nuclideCache, majorantCache, zaidCache)

  ! Public procedures
  public :: init
  public :: kill


contains


  !!
  !! Initialise Cache to given space
  !!
  !! Args:
  !!   Nmat [in]  -> Number of materials
  !!   Nnuc [in]  -> Number of nuclides
  !!   Nmaj [in]  -> Optional. Number of majorant Xss (Default = 1)
  !!   Nzaid [in] -> Optional. Number of nuclides with same ZAID
  !!
  !! Errors:
  !!   fatalError if Nmat, Nnuc or Nmaj is not a +ve value
  !!
  subroutine init(Nmat, Nnuc, Nmaj, Nzaid)
    integer(shortInt), intent(in)          :: Nmat
    integer(shortInt), intent(in)          :: Nnuc
    integer(shortInt),optional, intent(in) :: Nmaj
    integer(shortInt),optional, intent(in) :: Nzaid
    integer(shortInt)                      :: Nloc
    character(100),parameter :: Here = 'init (ceNeutronCache_mod.f90)'

    ! Make sure memory is clean
    call kill()

    ! Read default value of majorant XSs
    if(present(Nmaj)) then
      Nloc = Nmaj
    else
      Nloc = 1
    end if

    ! Chack the provided data
    if(Nmat < 1) call fatalError(Here,'Number of materials must be +ve! Not: '//numToChar(Nmat))
    if(Nnuc < 1) call fatalError(Here,'Number of nuclides must be +ve! Not: '//numToChar(Nmat))
    if(Nloc < 1) call fatalError(Here,'Number of majorant XSs must be +ve! Not: '//numToChar(Nmat))

    ! Allocate space
    ! Need to do in parallel region to allocate each copy
    !$omp parallel
    allocate(materialCache(Nmat))
    allocate(nuclideCache(Nnuc))
    allocate(majorantCache(Nloc))

    if(present(Nzaid)) then
      if (Nzaid > 0) then
        allocate(zaidCache(Nzaid))
      else
        call fatalError(Here,'Number of zaids must be +ve! Not: '//numToChar(Nzaid))
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
    if(allocated(materialCache)) deallocate (materialCache)
    if(allocated(nuclideCache))  deallocate (nuclideCache)
    if(allocated(majorantCache)) deallocate (majorantCache)
    if(allocated(zaidCache)) deallocate (zaidCache)
    !$omp end parallel
  end subroutine kill


end module ceNeutronCache_mod
