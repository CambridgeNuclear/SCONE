module materialSource_iTest

  use dictionary_class,     only : dictionary
  use dictParser_func,      only : charToDict
  use funit
  use geometry_inter,       only : geometry
  use geometryFactory_func, only : new_geometry
  use geometryReg_mod,      only : gr_geomIdx => geomIdx, gr_geomPtr => geomPtr, &
                                   gr_kill => kill
  use materialMenu_mod,     only : mm_kill => kill, mm_matIdx => matIdx
  use materialSource_class, only : materialSource
  use nuclearDataReg_mod,   only : ndReg_init => init, ndReg_activate => activate, &
                                   ndReg_kill => kill
  use numPrecision
  use particle_class,       only : particleState
  use RNG_class,            only : RNG
  use universalVariables,   only : P_NEUTRON_MG

  implicit none

  ! Geometry: a 4 x 4 x 4 cm box of mat1 with a void sphere of radius 1 cm at the
  ! centre. The source samples over the full geometry bounds, so the bounding box
  ! contains void (~6.5 % of its volume).
  character(*), parameter :: GEOMETRY_DEF = "                                    &
  &type geometryStd;                                                             &
  &boundary (0 0 0 0 0 0);                                                       &
  &graph { type shrunk; }                                                        &
  &surfaces {                                                                    &
  &  bound { id 1; type box; origin (0.0 0.0 0.0); halfwidth (2.0 2.0 2.0); }    &
  &  hole { id 2; type sphere; origin (0.0 0.0 0.0); radius 1.0; }               &
  &}                                                                             &
  &cells {                                                                       &
  &  void { id 1; type simpleCell; surfaces (-2); filltype mat; material void; } &
  &  mat { id 2; type simpleCell; surfaces (2); filltype mat; material mat1; }   &
  &}                                                                             &
  &universes {                                                                   &
  &  root { id 1; type rootUniverse; border 1; fill u<2>; }                      &
  &  uni { id 2; type cellUniverse; cells (1 2); }                               &
  &}"

  ! Nuclear data: a single MG material (4-group fixture shared with
  ! baseMgNeutronDatabase_iTest).
  character(*), parameter :: NUCLEAR_DATA_DEF = "                 &
  &handles { mg { type baseMgNeutronDatabase; PN P0; } }          &
  &materials {                                                    &
  &  mat1 { temp 273;                                             &
  &         composition { 1001.03 5.028E-02; 8016.03 2.505E-02; } &
  &         xsFile ./IntegrationTestFiles/mgMat1;                 &
  &       }                                                       &
  &}"

  ! Source: sample in mat1. No boundingBox entry, so defaults to the full
  ! geometry bounds which include the central void sphere.
  character(*), parameter :: SOURCE_DEF = &
  & "source { type materialSource; data mg; mat mat1; G 1; }"

contains
  !!
  !! Regression test: materialSource must resample positions that fall in
  !! void instead of aborting.
  !!
  !! With the void rejection missing, sampleParticle passes VOID_MAT to
  !! getMaterial and aborts with 'Nuclear data did not return neutron
  !! material' as soon as the void sphere is hit, which for this fixture
  !! happens with near certainty within the first few samples.
  !!
@Test
  subroutine testMaterialSourceWithVoidInBoundingBox()
    character(nameLen)           :: geomName, handleName, matName
    class(geometry), pointer     :: geom
    integer(shortInt)            :: i, matIdx
    type(dictionary)             :: geometryDict, nuclearDataDict, sourceDict
    type(materialSource)         :: source
    type(particleState)          :: state
    type(RNG)                    :: rand
    integer(longInt), parameter  :: SEED = 1_longInt
    integer(shortInt), parameter :: N_SAMPLES = 1000
    real(defReal), parameter     :: R_VOID_TOL = ONE - 1.0E-9_defReal

    ! Build nuclear data and material menu.
    call charToDict(nuclearDataDict, NUCLEAR_DATA_DEF)
    call ndReg_init(nuclearDataDict)

    ! Build geometry (materialMenu must already be loaded).
    call charToDict(geometryDict, GEOMETRY_DEF)
    geomName = 'materialSourceGeom'
    call new_geometry(geometryDict, geomName, silent = .true.)
    geom => gr_geomPtr(gr_geomIdx(geomName))

    ! Activate MG data.
    handleName = 'mg'
    call ndReg_activate(P_NEUTRON_MG, handleName, geom % activeMats(), silent = .true.)

    ! Build source.
    call charToDict(sourceDict, SOURCE_DEF)
    call source % init(sourceDict % getDictPtr('source'), geom)

    ! Sample sites with a fixed seed.
    call rand % init(SEED)
    matName = 'mat1'
    matIdx = mm_matIdx(matName)

    do i = 1, N_SAMPLES
      state = source % sampleParticle(rand)

      ! Every site must be in the requested material.
      @assertEqual(matIdx, state % matIdx)

      ! Every site must be outside the void sphere.
      @assertTrue(R_VOID_TOL < norm2(state % r), 'Source site inside the void sphere.')

    end do

    ! Clean up.
    call source % kill()
    call ndReg_kill()
    call gr_kill()
    call mm_kill()
    call geometryDict % kill()
    call nuclearDataDict % kill()
    call sourceDict % kill()

  end subroutine testMaterialSourceWithVoidInBoundingBox

end module materialSource_iTest