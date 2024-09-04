module urrProbabilityTables_iTest

  use numPrecision
  use dictionary_class,         only : dictionary
  use dictParser_func,          only : charToDict
  use particle_class,           only : particle
  use aceNeutronDatabase_class, only : aceNeutronDatabase
  use nuclearDatabase_inter,    only : nuclearDatabase
  use ceNeutronNuclide_inter,   only : ceNeutronNuclide, ceNeutronNuclide_CptrCast
  use aceNeutronNuclide_class,  only : aceNeutronNuclide, aceNeutronNuclide_CptrCast
  use neutronXSPackages_class,  only : neutronMicroXSs
  use materialMenu_mod,         only : mm_init => init
  use ceNeutronCache_mod,       only : zaidCache, nuclideCache
  use funit

  implicit none

  ! Material definitions
  character(*),parameter :: MAT_INPUT_STR = &
  & " uo2  {                   &
  &        composition {       &
  &        92235.03 1.0E-3;    &
  &        8016.03  2.0E-3;    &
  &        }                   &
  &      }"

  ! CE Neutron Database specification
  character(*),parameter :: ACE_INPUT_STR = &
  & "aceLibrary ./IntegrationTestFiles/testLib; ures 1 ; majorant 1; "

contains

  !!
  !! Test the use of probability tables
  !!
@Test
  subroutine test_urrProbabilityTables()
    type(aceNeutronDatabase), target  :: data
    class(nuclearDatabase), pointer   :: ptr
    type(dictionary)                  :: matDict
    type(dictionary)                  :: dataDict
    class(aceNeutronNuclide), pointer :: U235, O16
    real(defReal), dimension(3)       :: val
    real(defReal), dimension(2)       :: eBounds
    class(ceNeutronNuclide), pointer  :: nuc
    type(particle)                    :: p
    type(neutronMicroXSs)             :: microXSs
    real(defReal), parameter          :: TOL = 1.0E-6

    ! Prepare dictionaries
    call charToDict(matDict, MAT_INPUT_STR)
    call charToDict(dataDict, ACE_INPUT_STR)

    ! Build material menu
    call mm_init(matDict)

    ! Initialise data
    ptr => data
    call data % init(dataDict, ptr, silent = .true.)
    call data % activate([1], silent = .true.)

    !!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    !! Perform tests
    !!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Get nuclides
    U235  => aceNeutronNuclide_CptrCast( data % getNuclide(1))
    O16   => aceNeutronNuclide_CptrCast( data % getNuclide(2))

    !<><><><><><><><><><><><><><><><><><><><><><><><>
    ! Test probability tables

    @assertTrue(U235 % hasProbTab)
    @assertFalse(O16 % hasProbTab)

    !<><><><><><><><><><><><><><><><><><><><>
    ! Test energy bounds

    eBounds = U235 % probTab % getEbounds()

    @assertEqual(2.25E-3_defReal, eBounds(1), TOL)
    @assertEqual(2.5E-2_defReal,  eBounds(2), TOL)

    @assertEqual(O16 % urrE(1), ZERO)
    @assertEqual(O16 % urrE(2), ZERO)

    !<><><><><><><><><><><><><><><><><><><><><><><><>
    ! Test sampling from tables

    call U235 % probTab % sampleXSs(9.1E-3_defReal, 0.347_defReal, val)

    @assertEqual(0.98499622_defReal, val(1), TOL)
    @assertEqual(0.83939802_defReal, val(2), TOL)
    @assertEqual(0.8515398_defReal,  val(3), TOL)

    !<><><><><><><><><><><><><><><><><><><><><><><><>
    ! Test getting XSs

    ! U-235
    nuc  => ceNeutronNuclide_CptrCast( data % getNuclide(1))
    zaidCache(1) % E  = 9.1E-3_defReal
    zaidCache(1) % xi = 0.347_defReal
    nuclideCache(1) % E_tot = ONE

    call nuc % getMicroXSs(microXSs, 9.1E-3_defReal, ZERO, p % pRNG)

    @assertEqual(ONE, 15.317184903738868_defReal/ microXSs % total,            TOL)
    @assertEqual(ONE, 11.662135262310867_defReal/ microXSs % elasticScatter,   TOL)
    @assertEqual(ONE, 0.5743300000E-5_defReal   / microXSs % inelasticScatter, TOL)
    @assertEqual(ONE, 0.999051523404001_defReal / microXSs % capture,          TOL)
    @assertEqual(ONE, 2.655992374724002_defReal / microXSs % fission,          TOL)
    @assertEqual(ONE, 6.462838469906821_defReal / microXSs % nuFission,        TOL)

  end subroutine test_urrProbabilityTables


end module urrProbabilityTables_iTest
