program omp_app
  use numPrecision
  use openmp_func
  use endfConstants
  use universalVariables
  use genericProcedures,        only : linFind
  use dictionary_class,         only : dictionary
  use dictParser_func,          only : charToDict
  use charMap_class,            only : charMap
  use particle_class,           only : particle
  use aceNeutronDatabase_class, only : aceNeutronDatabase
  use nuclearDatabase_inter,    only : nuclearDatabase
  use ceNeutronMaterial_class,  only : ceNeutronMaterial, ceNeutronMaterial_TptrCast
  use ceNeutronNuclide_inter,   only : ceNeutronNuclide, ceNeutronNuclide_CptrCast
  use aceNeutronNuclide_class,  only : aceNeutronNuclide, aceNeutronNuclide_TptrCast
  use neutronXSPackages_class,  only : neutronMicroXSs, neutronMacroXSs
  use materialMenu_mod,         only : mm_init => init, mm_kill => kill

  implicit none

  ! Material definitions
  character(*),parameter :: MAT_INPUT_STR = " &
  & water { temp 273;                           &
  &       composition {                        &
  &       1001.03 5.028E-02;                   &
  &       8016.03 2.505E-02;                   &
  &       }                                    &
  &     }                                      &
  &uo2  { temp 1;                              &
  &        composition {                       &
  &        92233.03 2.286E-02;                 &
  &        8016.03  4.572E-02;                 &
  &        }                                   &
  &      }"

  ! CE Neutron Database specification
  character(*),parameter :: ACE_INPUT_STR = "aceLibrary ./IntegrationTestFiles/testLib;    "
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  type(aceNeutronDatabase), target :: data
  class(nuclearDatabase), pointer  :: ptr
  type(dictionary)                 :: matDict
  type(dictionary)                 :: dataDict
  integer(shortInt)                :: i
  real(defReal)                    :: E_MIN = 1.0E-9_defReal, E_MAX = 20.0_defReal
  integer(shortInt)                :: N = 1000
  real(defReal)                    :: E
  type(particle)                   :: p

  ! Prepare dictionaries
  call charToDict(matDict, MAT_INPUT_STR)
  call charToDict(dataDict, ACE_INPUT_STR)

  ! Build material menu
  call mm_init(matDict)

  ! Initialise data
  ptr => data
  call data % init(dataDict, ptr, silent = .true.)
  call data % activate([1,2])

  ! Interpolate XSs to verify

  call ompSetNumThreads(4)
  !$omp parallel do private(E,i,p)
  do i = 1, N
    E = (i-1) * log(E_MIN/E_MAX) / N
    E = exp(E) * E_MAX

    call p % build([ZERO, ZERO, ZERO], [ONE, ZERO, ZERO], E, ONE)

    print *, E, data % getTotalMatXS(p , 2), ompGetThreadNum()

  end do
  !$omp end parallel do

  !
  ! ! Interpolate XSs using number of threads
  !
  ! call ompSetNumThreads(2)
  !
  ! !$omp parallel
  ! if(ompGetThreadNum() == 0) then
  !   print *, "Master thread. # of threads active ", ompGetNumThreads(), "of", ompGetMaxThreads()
  ! end if
  !
  ! print *, 'Hello from thread: ', ompGetThreadNum()
  ! !$omp end parallel

end program omp_app
