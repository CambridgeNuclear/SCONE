module geometryStd_iTest

  use numPrecision
  use universalVariables
  use dictionary_class,  only : dictionary
  use dictParser_func,   only : charToDict
  use charMap_class,     only : charMap
  use dictParser_func,   only : fileToDict
  use coord_class,       only : coordList
  use geometry_inter,    only : geometry
  use geometryStd_class, only : geometryStd
  use visualiser_class,  only : visualiser
  use pFUnit_mod

  implicit none


contains

  !!
  !! Geometry integration test -> Simple 2x2 lattice
  !!
@Test
  subroutine test_lattice_geom()
    class(geometryStd), target, allocatable   :: geom
    character(*), parameter     :: path = './IntegrationTestFiles/Geometry/test_lat'
    type(charMap)               :: mats
    integer(shortInt)           :: i, idx, matIdx, uniqueID, event
    type(dictionary)            :: dict
    real(defReal), dimension(3) :: r, u, r_ref, u_ref
    type(dictionary),pointer    :: tempDict
    character(nameLen)          :: name
    type(coordList)             :: coords
    character(nameLen), dimension(:), allocatable :: keys
    integer(shortInt), dimension(10,10)           :: img
    real(defReal), dimension(6)                   :: aabb
    real(defReal)                                 :: maxDist
    class(geometry), pointer                      :: geomP
    type(visualiser)                              :: viz
    type(dictionary)                              :: vizDict
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    ! Load dictionary
    call fileToDict(dict, path)

    ! Load materials
    tempDict => dict % getDictPtr('nuclearData')
    tempDict => tempDict % getDictPtr('materials')
    call tempDict % keys(keys, 'dict')
    do i = 1, size(keys)
      call mats % add(keys(i), i)
    end do

    ! Build geometry
    call geom % init(dict, mats, silent=.true.)

    ! Get material at few locations
    name = 'water'
    idx = mats % get(name)
    call geom % whatIsAt(matIdx, uniqueID, [ZERO, ZERO, ZERO])
    @assertEqual(idx, matIdx)

    name = 'mox43'
    idx = mats % get(name)
    r = [0.63_defReal, -0.09_defReal, 0.0_defReal]
    u = [ZERO, -ONE, ZERO]
    call geom % whatIsAt(matIdx, uniqueID, r, u)
    @assertEqual(idx, matIdx)

    ! Place coordinates
    r = [0.1_defReal, 0.1_defReal, 0.0_defReal]
    u = [ZERO, ZERO, ONE]
    call coords % init(r, u)
    call geom % placeCoord(coords)

    ! Verify positions
    @assertEqual(r, coords % lvl(1) % r, TOL)
    @assertEqual(r, coords % lvl(2) % r, TOL)
    @assertEqual(r - [0.63_defReal, 0.63_defReal, 0.0_defReal], coords % lvl(3) % r, TOL)

    ! Verify directions
    @assertEqual(u, coords % lvl(1) % dir, TOL)
    @assertEqual(u, coords % lvl(2) % dir, TOL)
    @assertEqual(u, coords % lvl(3) % dir, TOL)

    ! Construct visualiser and verify slice plotting
    geomP => geom
    call charToDict(vizDict, ' ')
    call viz % init(geomP, vizDict)

    ! Slice plot -> Material
    call viz % slicePlot(img, [ZERO, ZERO, ZERO], 'z', 'material')

    ! Verify some pixels
    name = 'water'
    idx = mats % get(name)
    @assertEqual(idx, img(1, 1))
    @assertEqual(idx, img(2, 6))

    name = 'mox43'
    idx = mats % get(name)
    @assertEqual(idx, img(3, 7))

    name = 'uox'
    idx = mats % get(name)
    @assertEqual(idx, img(3, 3))

    ! Slice plot -> UniqueID
    r = [-0.63_defReal, -0.63_defReal, 0.0_defReal]
    call viz % slicePlot(img, r, 'z', 'uniqueID', [1.26_defReal, 1.26_defReal])

    ! Verify some pixels
    ! Note that this test depends on universe leyout order in geomGraph
    ! If it changes this test fill fail
    @assertEqual(2, img(5,5))
    @assertEqual(3, img(1,1))

    ! Verify bounds
    aabb = [-1.26_defReal, -1.26_defReal, 0.0_defReal, 1.26_defReal, 1.26_defReal, 0.0_defReal]
    @assertEqual(aabb, geom % bounds(), TOL)

    !*** Test teleport movement
    r = [ZERO, ZERO, ZERO]
    u = [-ONE, -TWO, ZERO]
    u = u/norm2(u)
    call coords % init(r, u)

    call geom % teleport(coords, 3.0_defReal)

    r_ref = [-1.1783592_defReal, -0.1632816_defReal, ZERO]
    u_ref = [ONE, -TWO, ZERO]
    u_ref = u_ref / norm2(u_ref)
    name = 'water'
    idx = mats % get(name)

    @assertEqual(r_ref, coords % lvl(1) % r, TOL)
    @assertEqual(u_ref, coords % lvl(1) % dir, TOL)
    @assertEqual(idx, coords % matIdx)

    !*** Test global movement
    r = [ZERO, ZERO, ZERO]
    u = [ZERO, -ONE, ZERO]
    call coords % init(r, u)

    ! Collosion movement
    maxDist = 1.0_defReal
    call geom % moveGlobal(coords, maxDist, event)

    r_ref = [ZERO, -1.0_defReal, ZERO]
    u_ref = u
    name = 'water'
    idx = mats % get(name)

    @assertEqual(r_ref, coords % lvl(1) % r, TOL)
    @assertEqual(u_ref, coords % lvl(1) % dir, TOL)
    @assertEqual(COLL_EV, event)
    @assertEqual(idx, coords % matIdx)
    @assertEqual(1.0_defReal, maxDist, TOL)

    ! Boundary Hit
    maxDist = 1.0_defReal
    call geom % moveGlobal(coords, maxDist, event)

    r_ref = [ZERO, 1.26_defReal, ZERO]
    u_ref = u_ref
    name = 'water'
    idx = mats % get(name)

    @assertEqual(r_ref, coords % lvl(1) % r, TOL)
    @assertEqual(u_ref, coords % lvl(1) % dir, TOL)
    @assertEqual(BOUNDARY_EV, event)
    @assertEqual(idx, coords % matIdx)
    @assertEqual(0.26_defReal, maxDist, TOL)

    !*** Normal Movment (easy case)
    r = [-0.63_defReal, -0.63_defReal, 0.0_defReal]
    u = [ZERO, -ONE, ZERO]
    call coords % init(r, u)
    call geom % placeCoord(coords)

    ! Local cell crossing
    maxDist = 1.0_defReal
    call geom % move(coords, maxDist, event)

    r_ref = [-0.63_defReal, -1.13_defReal, ZERO]
    u_ref = u_ref
    name = 'water'
    idx = mats % get(name)

    @assertEqual(r_ref, coords % lvl(1) % r, TOL)
    @assertEqual(u_ref, coords % lvl(1) % dir, TOL)
    @assertEqual(CROSS_EV, event)
    @assertEqual(idx, coords % matIdx)
    @assertEqual(0.5_defReal, maxDist, TOL)

    ! Boundary Hit
    maxDist = 1.0_defReal
    call geom % move(coords, maxDist, event)

    r_ref = [-0.63_defReal, 1.26_defReal, ZERO]
    u_ref = u_ref
    name = 'water'
    idx = mats % get(name)

    @assertEqual(r_ref, coords % lvl(1) % r, TOL)
    @assertEqual(u_ref, coords % lvl(1) % dir, TOL)
    @assertEqual(BOUNDARY_EV, event)
    @assertEqual(idx, coords % matIdx)
    @assertEqual(0.13_defReal, maxDist, TOL)

    ! Collision
    maxDist = 0.08_defReal
    call geom % move(coords, maxDist, event)

    r_ref = [-0.63_defReal, 1.18_defReal, ZERO]
    u_ref = u_ref
    name = 'water'
    idx = mats % get(name)

    @assertEqual(r_ref, coords % lvl(1) % r, TOL)
    @assertEqual(u_ref, coords % lvl(1) % dir, TOL)
    @assertEqual(COLL_EV, event)
    @assertEqual(idx, coords % matIdx)
    @assertEqual(0.08_defReal, maxDist, TOL)

    ! Kill geometry
    call geom % kill()

  end subroutine test_lattice_geom

  !!
  !! Test geometry with tilted cylinder
  !!
@Test
  subroutine test_tilted_cylinder()
    class(geometryStd), target, allocatable   :: geom
    character(*), parameter     :: path = './IntegrationTestFiles/Geometry/test_cyl'
    type(charMap)               :: mats
    integer(shortInt)           :: idxW, idxF, i
    type(dictionary)            :: dict
    type(dictionary),pointer    :: tempDict
    character(nameLen)          :: name
    character(nameLen), dimension(:), allocatable :: keys
    integer(shortInt), dimension(20,20)    :: img
    integer(shortInt), dimension(20,20,20) :: img3
    class(geometry), pointer    :: geomP
    type(visualiser)            :: viz
    type(dictionary)            :: vizDict
    real(defReal), dimension(3) :: r

    ! Load dictionary
    call fileToDict(dict, path)

    ! Load materials
    tempDict => dict % getDictPtr('nuclearData')
    tempDict => tempDict % getDictPtr('materials')
    call tempDict % keys(keys, 'dict')
    do i = 1, size(keys)
      call mats % add(keys(i), i)
    end do

    ! Build geometry
    call geom % init(dict, mats, silent=.true.)

    ! Get fuel and water index
    name = 'water'
    idxW = mats % get(name)

    name = 'mox43'
    idxF = mats % get(name)
    
    ! Construct visualiser and verify slice plotting
    geomP => geom
    call charToDict(vizDict, ' ')
    call viz % init(geomP, vizDict)

    !*** Test slice normal to x & y
    ! X-axis at 1.0
    r = [1.0_defReal, 0.0_defReal, 0.0_defReal]
    call viz % slicePlot(img, r, 'x', 'material')

    ! Test some pixels
    @assertEqual(idxW, img(8, 11))
    @assertEqual(idxW, img(17, 3))
    @assertEqual(idxF, img(10, 10))
    @assertEqual(idxF, img(18, 1))
    
    ! Y-axis at 3.0
    r = [0.0_defReal, 3.0_defReal, 0.0_defReal]
    call viz % slicePlot(img, r, 'y', 'material')

    @assertEqual(idxW, img(15, 1))
    @assertEqual(idxW, img(13, 4))
    @assertEqual(idxF, img(13, 3))
    @assertEqual(idxF, img(14, 2))

    !*** Test voxel plot
    ! Full plot
    ! Value of r is irrelevant
    call viz % voxelPlot(img3, r, 'material')

    ! Checksome against 2D plot
    r = [0.0_defReal, 2.75_defReal, 0.0_defReal]
    call viz % slicePlot(img, r, 'y', 'material')

    @assertEqual(img, img3(:,16,:))

    ! Small box all inside fuel
    r = [ 1.0_defReal, 0.0_defReal, 0.0_defReal]
    call viz % voxelPlot(img3, r, 'material', [0.5_defReal, 0.5_defReal, 0.5_defReal])

    @assertEqual(idxF, img3)

  end subroutine test_tilted_cylinder


end module geometryStd_iTest
