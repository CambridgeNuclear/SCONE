!!
!! Global Cache for MG Neutron Nuclear Data
!!
!! Any active nuclear database writes XSs to this module.
!! Any material handle reads from this module
!!
!! Public Members:
!!   matCache      -> Array of cached data for materials
!!   majorantCache -> Array of cached data for majorant XS
!!
!! NOTE:
!!   Cache arrays are deliberatly not targets. This is becouse there should be no pointers to the
!!   cache. Any update call can change energy of any value so it would not be possible that the
!!   energy of XSs pointed by pointers would not change silently.
!!
module mgNeutronCache_mod

  use numPrecision
  use genericProcedures,       only : fatalError, numToChar
  use neutronXsPackages_class, only : neutronMacroXSs
  use materialHandle_inter,    only : materialHandle

  implicit none
  private

  !!
  !! Structure that contains cached data for each MG Neutron Material
  !!
  !! Public Members:
  !!   G_tot  -> Energy group of the total XS in xss
  !!   G_tail -> Energy group of all XSs in xss except total
  !!   xss    -> Cached Cross-Section values
  !!
  type, public :: cacheMatDat
    real(defReal)         :: G_tot  = ZERO
    real(defReal)         :: G_tail = ZERO
    type(neutronMacroXSs) :: xss
    class(materialHandle), pointer :: mat
  end type cacheMatDat

  ! MEMBERS OF THE MODULE ARE GIVEN HERE
  type(cacheMatDat), dimension(:), allocatable, public   :: materialCache
  !$omp threadprivate(materialCache)

  ! Public procedures
  public :: init
  public :: kill


contains


  !!
  !! Initialise Cache to given space
  !!
  !! Args:
  !!   Nmat [in]  -> Number of materials
  !!   Nmaj [in]  -> Optional. Number of majorant Xss (Default = 1)
  !!
  !! Errors:
  !!   fatalError if Nmat or Nmaj is not a +ve value
  !!
  subroutine init(Nmat)
    integer(shortInt), intent(in)          :: Nmat
    character(100),parameter :: Here = 'init (ceNeutronCache_mod.f90)'

    ! Make sure memory is clean
    call kill()

    ! Chack the provided data
    if(Nmat < 1) call fatalError(Here,'Number of materials must be +ve! Not: '//numToChar(Nmat))

    ! Allocate space
    ! Need to do in parallel region to allocate each copy
    !$omp parallel
    allocate(materialCache(Nmat))
    !$omp end parallel

  end subroutine init

  !!
  !! Return Cache Module (Singleton) to uninitialised state
  !!
  subroutine kill()
    ! Need to deallocate on all threads
    !$omp parallel
    if(allocated(materialCache)) deallocate (materialCache)
    !$omp end parallel
  end subroutine kill


end module mgNeutronCache_mod
