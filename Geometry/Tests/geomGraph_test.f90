module geomGraph_test

  use numPrecision
  use universalVariables, only : OUTSIDE_MAT
  use dictionary_class,   only : dictionary
  use intMap_class,       only : intMap
  use uniFills_class,     only : uniFills
  use geomGraph_class,    only : geomGraph
  use funit

  implicit none


  ! Variables
  type(uniFills) :: geom

contains

  !!
  !! Set up fills
  !!
  !! Filling structure (->) contained materials (|) nested universes.
  !! See definitions for order in terms of localIDs
  !!
  !!  3 -> 0
  !!  | 7
  !!  | | 1001 -> 1 4
  !!  | | 1001 -> 1 4
  !!  | | 1003 -> 1 2
  !!
@Before
  subroutine set_up()
    integer(shortInt), dimension(:), allocatable :: fill
    type(intMap)                                 :: map

    call geom % init(7)

    ! Root uni id 3
    fill = [-7, OUTSIDE_MAT]
    call geom % addUniverse(1, 3, fill)
    call map % add(3, 1)

    ! Uni id 7
    fill = [-1001, -1001, -1003]
    call geom % addUniverse(2, 2, fill)
    call map % add(7, 2)

    ! Uni id 1001
    fill = [1, 4]
    call geom % addUniverse(3, 1001, fill)
    call map % add(1001, 3)

    ! Uni id 1002
    fill = [1, 3]
    call geom % addUniverse(4, 1002, fill)
    call map % add(1002, 4)

    ! Uni id 1003
    fill = [1, 2]
    call geom % addUniverse(5, 1003, fill)
    call map % add(1003, 5)

    ! Uni id 200
    fill = [-1001, -200]
    call geom % addUniverse(6, 200, fill)
    call map % add(200, 6)

    ! Uni id 201
    fill = [-1001, -200]
    call geom % addUniverse(7, 201, fill)
    call map % add(201, 7)

    ! Finish build
    call geom % finishBuild(map)
    call geom % setRoot(1)
  end subroutine set_up

  !!
  !! Clean up test enviroment
  !!
@After
  subroutine clean_up()

    call geom % kill()

  end subroutine clean_up

  !!
  !! Test geometry graphs build in 'shrunk' version
  !!
@Test
  subroutine test_shrunk()
    type(geomGraph)   :: graph
    integer(shortInt) :: idx, id
    type(dictionary)  :: dict

    ! Create input dictionary
    call dict % init(1)
    call dict % store('type', 'shrunk')

    call graph % kill() ! Should work fine for uninitialised as well
    call graph % init(geom, dict)

    ! Verify location array
    ! Test will fail if traverse order is change (even if structure is correct)!
    @assertEqual([-2, 0, -3, -3, -5, 1, 4, 1, 2], graph % array % idx)
    @assertEqual([ 3, 0,  6,  6,  8, 1, 2, 3, 4], graph % array % id)

    ! Verify number of unique cells
    @assertEqual(4, graph % uniqueCells)

    ! Verify used materials
    @assertEqual([1, 2, 4], graph % usedMats)

    ! Test gtting content
    ! Universe
    call graph % getFill(idx, id, 1, 1)
    @assertEqual(-2, idx)
    @assertEqual(3, id)

    call graph % getFill(idx, id, 3, 3)
    @assertEqual(-5, idx)
    @assertEqual(8, id)

    ! Material
    call graph % getFill(idx, id, 8, 2)
    @assertEqual(2, idx)
    @assertEqual(4, id)

    ! Clean up
    call graph % kill()

  end subroutine test_shrunk

  !!
  !! Test extended setting with instances copying
  !!
@Test
  subroutine test_extended()
    type(geomGraph)   :: graph
    integer(shortInt) :: idx, id
    type(dictionary)  :: dict

    ! Create input dictionary
    call dict % init(1)
    call dict % store('type', 'extended')

    call graph % kill() ! Should work fine for uninitialised as well
    call graph % init(geom, dict)

    ! Verify location array
    ! Test will fail if traverse order is change (even if structure is correct)!
    @assertEqual([-2, 0, -3, -3, -5, 1, 4, 1, 4, 1, 2], graph % array % idx)
    @assertEqual([ 3, 0,  6,  8, 10, 1, 2, 3, 4, 5, 6], graph % array % id)

    ! Verify number of unique cells
    @assertEqual(6, graph % uniqueCells)

    ! Verify used materials
    @assertEqual([1, 2, 4], graph % usedMats)

    ! Test gtting content
    ! Universe
    call graph % getFill(idx, id, 1, 1)
    @assertEqual(-2, idx)
    @assertEqual(3, id)

    call graph % getFill(idx, id, 3, 3)
    @assertEqual(-5, idx)
    @assertEqual(10, id)

    ! Material
    call graph % getFill(idx, id, 8, 2)
    @assertEqual(4, idx)
    @assertEqual(4, id)

    ! Clean up
    call graph % kill()

  end subroutine test_extended


end module geomGraph_test
