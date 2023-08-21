module thermalScatteringData_iTest

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
  use ceNeutronCache_mod,       only : nuclideCache
  use funit

  implicit none

  ! Material definitions
  character(*),parameter :: MAT_INPUT_STR =        &
  & "water { temp 1;                               &
  &       moder {1001.03 h-h2o.49; }               &
  &       composition {                            &
  &       1001.03  2.0E-3;                         &
  &       8016.03  1.0E-3;                         &
  &                   }                            &
  &      }                                         &
  &  graphite { temp 1;                                &
  &          moder {6012.06  grph30.46;}               &
  &          composition {                             &
  &          6012.06 2.0E-3;                           &
  &                       }                            &
  &            }"

  ! CE Neutron Database specification
  character(*),parameter :: ACE_INPUT_STR = &
  & "aceLibrary ./IntegrationTestFiles/testLib; "

contains

  !!
  !! Test the use of thermal scattering libraries
  !!
@Test
  subroutine test_thermalScatteringData()
    type(aceNeutronDatabase), target  :: data
    class(nuclearDatabase), pointer   :: ptr
    type(dictionary)                  :: matDict
    type(dictionary)                  :: dataDict
    class(aceNeutronNuclide), pointer :: H1, O16, C12
    real(defReal)                     :: val
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
    call data % activate(([1,2]))

    !!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    !! Perform tests
    !!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Get nuclides
    O16 => aceNeutronNuclide_CptrCast( data % getNuclide(1))
    H1  => aceNeutronNuclide_CptrCast( data % getNuclide(2))
    C12 => aceNeutronNuclide_CptrCast( data % getNuclide(3))

    !<><><><><><><><><><><><><><><><><><><><><><><><>
    ! Test scattering tables

    @assertTrue(H1 % hasThData)
    @assertFalse(O16 % hasThData)

    @assertFalse(H1 % thData % hasElastic)

    @assertTrue(C12 % hasThData)
    @assertTrue(C12 % thData % hasElastic)
    @assertTrue(C12 % thData % isCoherent)

    !<><><><><><><><><><><><><><><><><><><><>
    ! Test energy bounds
    eBounds = H1 % thData % getEbounds('inelastic')

    @assertEqual(1.000E-11_defReal, eBounds(1), TOL)
    @assertEqual(1.000E-5_defReal,  eBounds(2), TOL)

    @assertEqual(O16 % SabInel(1), ZERO)
    @assertEqual(O16 % SabInel(2), ZERO)

    eBounds = C12 % thData % getEbounds('elastic')

    @assertEqual(1.000E-11_defReal, eBounds(1), TOL)
    @assertEqual(4.9000E-06,  eBounds(2), TOL)

    !<><><><><><><><><><><><><><><><><><><><><><><><>
    ! Test sampling from tables

    val = H1 % thData % getInelXS(1.8E-6_defReal)
    @assertEqual(21.018654322_defReal, val, TOL)

    val = H1 % thData % getElXS(1.8E-6_defReal)
    @assertEqual(ZERO, val, TOL)

    !<><><><><><><><><><><><><><><><><><><><><><><><>
    ! Test Getting material XSs
    ! water

    call p % build([ZERO, ZERO, ZERO], [ONE, ZERO, ZERO], 1.0E-6_defReal, ONE)

    ! Total XS of water
    p % E = 1.8E-6_defReal
    @assertEqual(ONE, data % getTotalMatXS(p , 1) / 0.0459700882_defReal , TOL)

    !<><><><><><><><><><><><><><><><><><><><><><><><>
    ! Test getting XSs
    ! H-1
    nuc  => ceNeutronNuclide_CptrCast( data % getNuclide(2))
    nuclideCache(2) % E_tot = ONE
    nuclideCache(2) % needsSabInel = .true.

    call nuc % getMicroXSs(microXSs, 1.8E-6_defReal, p % pRNG)

    @assertEqual(ONE, 21.05810233858_defReal/ microXSs % total,          TOL)
    @assertEqual(ONE, 21.01865432_defReal / microXSs % inelasticScatter, TOL)
    @assertEqual(ONE, 3.94480160E-002_defReal / microXSs % capture,      TOL)
    @assertEqual(ZERO, microXSs % elasticScatter,   TOL)
    @assertEqual(ZERO, microXSs % fission,          TOL)
    @assertEqual(ZERO, microXSs % nuFission,        TOL)

  end subroutine test_thermalScatteringData


end module thermalScatteringData_iTest
