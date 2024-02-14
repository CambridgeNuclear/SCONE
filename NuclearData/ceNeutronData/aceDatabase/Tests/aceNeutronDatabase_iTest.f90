module aceNeutronDatabase_iTest

  use numPrecision
  use endfConstants
  use universalVariables
  use genericProcedures,        only : linFind
  use dictionary_class,         only : dictionary
  use dictParser_func,          only : charToDict
  use charMap_class,            only : charMap
  use particle_class,           only : particle
  use aceNeutronDatabase_class, only : aceNeutronDatabase
  use nuclearDatabase_inter,    only : nuclearDatabase
  use materialHandle_inter,     only : materialHandle
  use nuclideHandle_inter,      only : nuclideHandle
  use reactionHandle_inter,     only : reactionHandle
  use ceNeutronMaterial_class,  only : ceNeutronMaterial, ceNeutronMaterial_TptrCast
  use ceNeutronNuclide_inter,   only : ceNeutronNuclide, ceNeutronNuclide_CptrCast
  use aceNeutronNuclide_class,  only : aceNeutronNuclide, aceNeutronNuclide_TptrCast
  use neutronXSPackages_class,  only : neutronMicroXSs, neutronMacroXSs
  use materialMenu_mod,         only : mm_init => init, mm_kill => kill
  use pFUnit_mod

  implicit none

  ! Material definitions
  character(*),parameter :: MAT_INPUT_STR = &
  & "water { temp 273;           &
  &       composition {          &
  &       1001.03 5.028E-02;     &
  &       8016.03 2.505E-02;     &
  &                   }          &
  &        }                     &
  &  uo2  { temp 1;              &
  &        composition {         &
  &        92233.03 2.286E-02;   &
  &        8016.03  4.572E-02;   &
  &                    }         &
  &       }"

  ! CE Neutron Database specification
  character(*),parameter :: ACE_INPUT_STR = &
  & "aceLibrary ./IntegrationTestFiles/testLib; "

contains

  !!
  !! One big monster test to avoid expensive set up each test
  !!
@Test
  subroutine test_aceNeutronDatabase()
    type(aceNeutronDatabase), target :: data
    class(nuclearDatabase), pointer  :: ptr
    type(dictionary)                 :: matDict
    type(dictionary)                 :: dataDict
    type(ceNeutronMaterial),pointer  :: mat
    class(ceNeutronNuclide), pointer :: nuc
    type(aceNeutronNuclide), pointer :: nuc2
    class(reactionHandle), pointer   :: reac
    type(charMap), pointer           :: matNames
    character(nameLen)               :: name
    type(particle)                   :: p
    type(neutronMicroXSs)            :: microXSs
    type(neutronMacroXSs)            :: macroXSs
    real(defReal)                    :: t1, t2
    integer(shortInt)                :: i, H1, O16, U233
    real(defReal), parameter         :: TOL = 1.0E-6

    ! Prepare dictionaries
    call charToDict(matDict, MAT_INPUT_STR)
    call charToDict(dataDict, ACE_INPUT_STR)

    ! Build material menu
    call mm_init(matDict)

    ! Initialise data
    ptr => data
    call data % init(dataDict, ptr, silent = .true.)
    call data % activate([1,2], silent = .true.)

    !!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    !! Perform tests
    !!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    !!<><><><><><><><><><><><><><><><><><><><><><
    !! Test getting material
    mat => null()

    ! Get invalid materials
    @assertNotAssociated( ceNeutronMaterial_TptrCast( data % getMaterial(0)))
    @assertNotAssociated( ceNeutronMaterial_TptrCast( data % getMaterial(-4)))
    @assertNotAssociated( ceNeutronMaterial_TptrCast( data % getMaterial(3)))

    ! Get water
    mat => ceNeutronMaterial_TptrCast( data % getMaterial(1))
    @assertAssociated(mat)

    ! Make sure densities are present
    @assertTrue(targetNotFound /= linFind(mat % dens, 5.028E-02_defReal, TOL), "H-1 dens is absent")
    @assertTrue(targetNotFound /= linFind(mat % dens, 2.505E-02_defReal, TOL), "O-16 dens is absent")

    @assertEqual(1, mat % matIdx)
    @assertFalse( mat % isFissile())
    @assertAssociated(mat % data, data)

    ! Get UO2
    mat => ceNeutronMaterial_TptrCast( data % getMaterial(2))
    @assertAssociated(mat)

    @assertTrue(targetNotFound /= linFind(mat % dens, 2.286E-02_defReal, TOL), "U-233 dens is absent")
    @assertTrue(targetNotFound /= linFind(mat % dens, 4.572E-02_defReal, TOL), "O-16 dens is absent")

    @assertEqual(2, mat % matIdx)
    @assertTrue( mat % isFissile())
    @assertAssociated(mat % data, data)

    !<><><><><><><><><><><><><><><><><><><><><><><><>
    ! Test getting nuclides
    nuc => null()

    ! Get invalid Nuclides
    @assertNotAssociated( ceNeutronNuclide_CptrCast( data % getNuclide(0)))
    @assertNotAssociated( ceNeutronNuclide_CptrCast( data % getNuclide(4)))
    @assertNotAssociated( ceNeutronNuclide_CptrCast( data % getNuclide(-3)))

    ! Get Valid Nuclides
    @assertAssociated( ceNeutronNuclide_CptrCast( data % getNuclide(1)))
    @assertAssociated( ceNeutronNuclide_CptrCast( data % getNuclide(2)))
    @assertAssociated( ceNeutronNuclide_CptrCast( data % getNuclide(3)))

    !<><><><><><><><><><><><><><><><><><><><><><><><>
    ! Test getting reactions
    reac => null()


    ! Nuclides can have diffrent indexes if Hashes do not work correctly
    ! Need to explicitly find which index correcponds to which nuclide
    ! Usually will be the following
    ! Nuclides 1 -> O-16
    !          2 -> U-233
    !          3 -> H-1
    do i=1,3
      nuc2 => aceNeutronNuclide_TptrCast( data % getNuclide(i))
      select case(trim(adjustl(nuc2 % ZAID)))
        case('1001.03c')
          H1 = i
        case('8016.03c')
          O16 = i
        case('92233.03c')
          U233 = i
      end select
    end do

    ! Get Invalid Reaction
    @assertNotAssociated( data % getReaction(N_Nl(3), H1))
    @assertNotAssociated( data % getReaction(macroEscatter, O16))
    @assertNotAssociated( data % getReaction(N_total, U233))
    @assertNotAssociated( data % getReaction(N_N_elastic, 4))
    @assertNotAssociated( data % getReaction(N_N_elastic, -2))
    @assertNotAssociated( data % getReaction(N_N_elastic, 0))

    ! Get Valid Reaction
    @assertAssociated( data % getReaction(N_fission, U233))
    @assertAssociated( data % getReaction(N_N_elastic, O16))
    @assertAssociated( data % getReaction(N_NL(3), U233))

    !<><><><><><><><><><><><><><><><><><><><><><><><>
    ! Get material names dictionary
    matNames => data % matNamesMap()
    name = 'water'
    @assertTrue( 0 /= matNames % getOrDefault(name, 0))

    name = 'uo2'
    @assertTrue( 0 /= matNames % getOrDefault(name, 0))


    !<><><><><><><><><><><><><><><><><><><><><><><><>
    ! Test getting nuclide XSs
     !
    ! SET UP as a regression test! Take values with the grain of salt!
    !

    ! H-1
    nuc  => ceNeutronNuclide_CptrCast( data % getNuclide(H1))
    @assertEqual(ONE, 20.765855864000002_defReal/ nuc % getTotalXS(1.1E-6_defReal, p % pRNG), TOL)

    call nuc % getMicroXSs(microXSs, 5.6E-3_defReal, p % pRNG)

    ! Absent XSs
    @assertEqual(ZERO, microXSs % fission)
    @assertEqual(ZERO, microXSs % nuFission)
    @assertEqual(ZERO, microXSs % inelasticScatter)

    ! Present XSs
    @assertEqual(ONE, 19.731020820000000_defReal   / microXSs % total,          TOL)
    @assertEqual(ONE, 19.730326000000000_defReal   / microXSs % elasticScatter, TOL)
    @assertEqual(ONE, 6.948036800000000e-04_defReal/ microXSs % capture,        TOL)

    !<><><><><><><><><><><><><><><><><><><><><><><><>
    ! Test Getting material XSs
    !

    call p % build([ZERO, ZERO, ZERO], [ONE, ZERO, ZERO], 1.0E-6_defReal, ONE)

    ! Total XS of water
    p % E = 1.1E-6_defReal
    @assertEqual(ONE, data % getTotalMatXS(p , 1)/1.1406745607419302_defReal , TOL)

    p % E = 19.9_defReal
    @assertEqual(ONE, data % getTrackingXS(p, 1, MATERIAL_XS)/6.539039844E-02_defReal , TOL)


    ! Total XS of UO2
    p % E = 1.1E-6_defReal
    @assertEqual(ONE, data % getTotalMatXS(p , 2)/4.4149556129495560_defReal , TOL)

    p % E = 19.9_defReal
    @assertEqual(ONE, data % getTrackingXS(p , 2, MATERIAL_XS)/0.21869599644_defReal , TOL)

    ! Majorant
    p % E = 1.1E-6_defReal
    @assertEqual(ONE, data % getMajorantXS(p) /4.4149556129495560_defReal , TOL)
    @assertEqual(ONE, data % getTrackingXS(p , 3, MAJORANT_XS) /4.4149556129495560_defReal , TOL)

    p % E = 19.9_defReal
    @assertEqual(ONE, data % getMajorantXS(p)/0.21869599644_defReal , TOL)
    @assertEqual(ONE, data % getTrackingXS(p , 3, MAJORANT_XS) /0.21869599644_defReal , TOL)

    ! Check that results are the same with on-the-fly majorant
    data % hasMajorant = .false.

    p % E = 1.1E-6_defReal
    @assertEqual(ONE, data % getMajorantXS(p) /4.4149556129495560_defReal , TOL)

    p % E = 19.9_defReal
    @assertEqual(ONE, data % getMajorantXS(p)/0.21869599644_defReal , TOL)

    !<><><><><><><><><><><><><><><><><><><><>
    ! Test getting Macroscopic XSs
    !
    ! Water
    mat => ceNeutronMaterial_TptrCast( data % getMaterial(1))
    call mat % getMacroXSs(macroXss, 3.6E-1_defReal, p % pRNG)

    ! Absent XSs
    @assertEqual(ZERO, macroXSs % fission)
    @assertEqual(ZERO, macroXSs % nuFission)

    @assertEqual(ONE, 0.466713100775700_defReal/ macroXSs     % total, TOL)
    @assertEqual(ONE, 0.466710902790000_defReal/ macroXSs     % elasticScatter, TOL)
    @assertEqual(ZERO, macroXSs % inelasticScatter, TOL)
    @assertEqual(ONE, 2.198066842597500e-06_defReal/ macroXSs % capture, TOL)

    ! Water with some inelastic collisions
    call mat % getMacroXSs(macroXss, 6.525_defReal, p % pRNG)

    @assertEqual(ONE, macroXSs % inelasticScatter/1.903667536E-04_defReal, TOL)


    !<><><><><><><><><><><><><><><><><><><><>
    ! Test getting energy bounds
    !
    call data % energyBounds(t1,t2)
    @assertEqual(1.0E-11_defReal, t1, TOL)
    @assertEqual(20.0_defReal,    t2, TOL)

    ! Clean everything
    call data % kill()
    call mm_kill()

  end subroutine test_aceNeutronDatabase

end module aceNeutronDatabase_iTest
