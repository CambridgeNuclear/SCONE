module baseMgNeutronDatabase_iTest

  use numPrecision
  use endfConstants
  use pFUnit_mod
  use universalVariables
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use particle_class,     only : particle

  ! Nuclear Data Objects & Interfaces
  use baseMgNeutronDatabase_class, only : baseMgNeutronDatabase, baseMgNeutronDatabase_CptrCast, &
                                          baseMgNeutronDatabase_TptrCast
  use baseMgNeutronMaterial_class, only : baseMgNeutronMaterial, baseMgNeutronMaterial_CptrCast, &
                                          baseMgNeutronMaterial_TptrCast
  use fissionMG_class,             only : fissionMG, fissionMG_TptrCast
  use multiScatterMG_class,        only : multiScatterMG, multiScatterMG_CptrCast, &
                                          multiScatterMG_TptrCast
  use multiScatterP1MG_class,      only : multiScatterP1MG, multiScatterP1MG_TptrCast
  use materialMenu_mod,            only : mm_init => init, mm_kill => kill
  use nuclearDatabase_inter,       only : nuclearDatabase
  use materialHandle_inter,        only : materialHandle
  use nuclideHandle_inter,         only : nuclideHandle
  use neutronXsPackages_class,     only : neutronMacroXSs
  use reactionHandle_inter,        only : reactionHandle



  implicit none

  ! Material definitions
  character(*),parameter :: MAT_INPUT_STR = "   &
  mat1 { temp 273;                              &
         composition {                          &
         1001.03 5.028E-02;                     &
         8016.03 2.505E-02;                     &
         }                                      &
         xsFile ./IntegrationTestFiles/mgMat1;  &
       }                                        &
  mat2  { temp 1;                               &
          composition {                         &
          92233.03 2.286E-02;                   &
          8016.03  4.572E-02;                   &
          }                                     &
          xsFile ./IntegrationTestFiles/mgMat2; &
        }"


contains

  !!
  !! Monster test to build and verify data in baseMgNeutronDatabase with P0 scattering
  !!
@Test
  subroutine testBaseMgNeutronDatabaseWithP0()
    type(baseMgNeutronDatabase), target  :: database
    class(nuclearDatabase), pointer      :: data_ptr
    type(dictionary)                     :: databaseDef
    type(dictionary)                     :: matMenuDict
    type(particle)                       :: p
    type(neutronMacroXSs)                :: xss
    type(baseMgNeutronMaterial),pointer  :: mat
    class(baseMgNeutronMaterial),pointer :: matClass
    class(reactionHandle), pointer       :: reac
    real(defReal),parameter :: TOL = 1.0E-6_defReal


    data_ptr => database

    ! Load materialMenu
    call charToDict(matMenuDict, MAT_INPUT_STR)
    call mm_init(matMenuDict )

    ! Build database
    call databaseDef % init(1)
    call databaseDef % store('PN','P0')
    call database % init(databaseDef, data_ptr, silent = .true.)
    call database % activate([1], silent = .true.)

    ! Varify number of groups
    @assertEqual(4, database % nGroups())

    ! Test getting Transport XS
    p % G = 1
    @assertEqual(2.1_defReal, database % getTrackingXS(p, 1, MATERIAL_XS), TOL)

    ! Test getting Total XS
    p % G = 1
    @assertEqual(3.1_defReal, database % getTotalMatXS(p, 2), TOL)

    p % G = 3
    @assertEqual(6.0_defReal, database % getTotalMatXS(p, 1), TOL)

    ! Test getting Majorant
    p % G = 1
    @assertEqual(2.1_defReal, database % getMajorantXS(p), TOL)
    @assertEqual(2.1_defReal, database % getTrackingXS(p, 1, MAJORANT_XS), TOL)

    ! Get a material and verify macroXSS
    mat => baseMgNeutronMaterial_TptrCast(database % getMaterial(2))
    @assertTrue(associated(mat), "Type Ptr Cast has failed")
    call mat % getMacroXSs(xss, 1, p % pRNG)

    ! Check that is fissile
    @assertTrue(mat % isFissile(), "Is not fissile but should")

    @assertEqual(3.1_defReal, xss % total, TOL)
    @assertEqual(ZERO, xss % elasticScatter, TOL)
    @assertEqual(1.1_defReal, xss % inelasticScatter, TOL)
    @assertEqual(1.0_defReal, xss % capture, TOL)
    @assertEqual(1.0_defReal, xss % fission, TOL)
    @assertEqual(2.3_defReal, xss % nuFission, TOL)

    matClass => baseMgNeutronMaterial_CptrCast(database % getMaterial(1))
    @assertTrue(associated(matClass), "Type Ptr Cast has failed")
    call matClass % getMacroXSs(xss, 4, p % pRNG)

    @assertFalse(matClass % isFissile(), "Is fissile but should not")

    @assertEqual(7.1_defReal, xss % total, TOL)
    @assertEqual(ZERO, xss % elasticScatter, TOL)
    @assertEqual(3.1_defReal, xss % inelasticScatter, TOL)
    @assertEqual(4.0_defReal, xss % capture, TOL)
    @assertEqual(0.0_defReal, xss % fission, TOL)
    @assertEqual(0.0_defReal, xss % nuFission, TOL)

    ! Get some invalid Materials
    mat => baseMgNeutronMaterial_TptrCast(database % getMaterial(0))
    @assertFalse(associated(mat))

    mat => baseMgNeutronMaterial_TptrCast(database % getMaterial(-2))
    @assertFalse(associated(mat))

    mat => baseMgNeutronMaterial_TptrCast(database % getMaterial(3))
    @assertFalse(associated(mat))

    ! Get Fission reaction and verify type
    reac => fissionMG_TptrCast(database % getReaction(macroFission, 1))
    @assertFalse(associated(reac), "Pointer for the mission reaction is not null")

    reac => fissionMG_TptrCast(database % getReaction(macroFission, 2))
    @assertTrue(associated(reac), "Pointer fission reaction is wrong type or null")

    ! Get Scattering reaction and verify type
    reac => multiScatterMG_TptrCast(database % getReaction(macroIEScatter, 1))
    @assertTrue(associated(reac), "Wrong type of scattering reaction")

    ! Get some invalid reactions
    reac => database % getReaction(anyScatter, 0)
    @assertFalse(associated(reac))

    reac => database % getReaction(anyScatter, -1)
    @assertFalse(associated(reac))

    reac => database % getReaction(anyScatter, 3)
    @assertFalse(associated(reac))

    reac => database % getReaction(anyCapture, 1)
    @assertFalse(associated(reac))

    ! **** Note that anyFission is not present !
    reac => database % getReaction(anyFission, 2)
    @assertFalse(associated(reac))

    ! Test getting nuclide
    @assertFalse(associated(database % getNuclide(1)))

    ! Clean up
    call database % kill()
    call mm_kill()
    call matMenuDict % kill()
    call databaseDef % kill()

  end subroutine testBaseMgNeutronDatabaseWithP0

  !!
  !! Monster test to build and verify data in baseMgNeutronDatabase with P1 scattering
  !! *Copy and pasted from the above with only the type of scattering changed
  !!
@Test
  subroutine testBaseMgNeutronDatabaseWithP1()
    type(baseMgNeutronDatabase), target  :: database
    class(nuclearDatabase), pointer      :: data_ptr
    type(dictionary)                     :: databaseDef
    type(dictionary)                     :: matMenuDict
    type(particle)                       :: p
    type(neutronMacroXSs)                :: xss
    type(baseMgNeutronMaterial),pointer  :: mat
    class(baseMgNeutronMaterial),pointer :: matClass
    class(reactionHandle), pointer       :: reac
    real(defReal),parameter :: TOL = 1.0E-6_defReal


    data_ptr => database

    ! Load materialMenu
    call charToDict(matMenuDict, MAT_INPUT_STR)
    call mm_init(matMenuDict )

    ! Build database
    call databaseDef % init(1)
    call databaseDef % store('PN','P1')
    call database % init(databaseDef, data_ptr, silent = .true.)
    call database % activate([1], silent = .true.)

    ! Varify number of groups
    @assertEqual(4, database % nGroups())

    ! Test getting Transport XS
    p % G = 1
    @assertEqual(2.1_defReal, database % getTrackingXS(p, 1, MATERIAL_XS), TOL)

    ! Test getting Total XS
    p % G = 1
    @assertEqual(3.1_defReal, database % getTotalMatXS(p, 2), TOL)

    p % G = 3
    @assertEqual(6.0_defReal, database % getTotalMatXS(p, 1), TOL)

    ! Test getting Majorant
    p % G = 1
    @assertEqual(2.1_defReal, database % getMajorantXS(p), TOL)
    @assertEqual(2.1_defReal, database % getTrackingXS(p, 1, MAJORANT_XS), TOL)

    ! Get a material and verify macroXSS
    mat => baseMgNeutronMaterial_TptrCast(database % getMaterial(2))
    @assertTrue(associated(mat), "Type Ptr Cast has failed")
    call mat % getMacroXSs(xss, 1, p % pRNG)

    ! Check that is fissile
    @assertTrue(mat % isFissile(), "Is not fissile but should")

    @assertEqual(3.1_defReal, xss % total, TOL)
    @assertEqual(ZERO, xss % elasticScatter, TOL)
    @assertEqual(1.1_defReal, xss % inelasticScatter, TOL)
    @assertEqual(1.0_defReal, xss % capture, TOL)
    @assertEqual(1.0_defReal, xss % fission, TOL)
    @assertEqual(2.3_defReal, xss % nuFission, TOL)

    matClass => baseMgNeutronMaterial_CptrCast(database % getMaterial(1))
    @assertTrue(associated(matClass), "Type Ptr Cast has failed")
    call matClass % getMacroXSs(xss, 4, p % pRNG)

    @assertFalse(matClass % isFissile(), "Is fissile but should not")

    @assertEqual(7.1_defReal, xss % total, TOL)
    @assertEqual(ZERO, xss % elasticScatter, TOL)
    @assertEqual(3.1_defReal, xss % inelasticScatter, TOL)
    @assertEqual(4.0_defReal, xss % capture, TOL)
    @assertEqual(0.0_defReal, xss % fission, TOL)
    @assertEqual(0.0_defReal, xss % nuFission, TOL)

    ! Get some invalid Materials
    mat => baseMgNeutronMaterial_TptrCast(database % getMaterial(0))
    @assertFalse(associated(mat))

    mat => baseMgNeutronMaterial_TptrCast(database % getMaterial(-2))
    @assertFalse(associated(mat))

    mat => baseMgNeutronMaterial_TptrCast(database % getMaterial(3))
    @assertFalse(associated(mat))

    ! Get Fission reaction and verify type
    reac => fissionMG_TptrCast(database % getReaction(macroFission, 1))
    @assertFalse(associated(reac), "Pointer for the mission reaction is not null")

    reac => fissionMG_TptrCast(database % getReaction(macroFission, 2))
    @assertTrue(associated(reac), "Pointer fission reaction is wrong type or null")

    ! Get Scattering reaction and verify type
    reac => multiScatterP1MG_TptrCast(database % getReaction(macroIEScatter, 1))
    @assertTrue(associated(reac), "Wrong type of scattering reaction")

    ! Get some invalid reactions
    reac => database % getReaction(anyScatter, 0)
    @assertFalse(associated(reac))

    reac => database % getReaction(anyScatter, -1)
    @assertFalse(associated(reac))

    reac => database % getReaction(anyScatter, 3)
    @assertFalse(associated(reac))

    reac => database % getReaction(anyCapture, 1)
    @assertFalse(associated(reac))

    ! **** Note that anyFission is not present !
    reac => database % getReaction(anyFission, 2)
    @assertFalse(associated(reac))

    ! Test getting nuclide
    @assertFalse(associated(database % getNuclide(1)))

    ! Clean up
    call database % kill()
    call mm_kill()
    call matMenuDict % kill()
    call databaseDef % kill()

  end subroutine testBaseMgNeutronDatabaseWithP1



end module baseMgNeutronDatabase_iTest
