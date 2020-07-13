module uniFills_test

  use numPrecision
  use genericProcedures,  only : linFind
  use universalVariables, only : OUTSIDE_MAT, targetNotFound
  use intMap_class,       only : intMap
  use uniFills_class,     only : uniFills
  use pFUnit_mod

  implicit none

  ! Variables
  type(uniFills) :: geom

contains

  !!
  !! Setup an acyclic geometry
  !!
  !! 7 universes
  !! 5 used
  !! 3 level nesting
  !! Contains outside below root
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
    fill = [-1001, -1002, -1003, -1001]
    call geom % addUniverse(2, 2, fill)
    call map % add(7, 2)

    ! Uni id 1001
    fill = [1, OUTSIDE_MAT]
    call geom % addUniverse(3, 1001, fill)
    call map % add(1001, 3)

    ! Uni id 1002
    fill = [1, 2]
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
  !! Clean up
  !!
@After
  subroutine clean_up()
    call geom % kill()

  end subroutine clean_up

  !!
  !! Test check for cycles
  !!
  !! Does not use the common graph
  !!
@Test
  subroutine test_cycles()
    type(uniFills)                               :: graph
    integer(shortInt), dimension(:), allocatable :: fill
    type(intMap)                                 :: map

    ! Build a geometry graph with recursion in the same universe (loop in 2nd)
    call graph % init(3)

    fill = [-2, OUTSIDE_MAT]
    call graph % addUniverse(1, 1, fill)
    call map % add(1, 1)

    fill = [-2, -3, 2, 3]
    call graph % addUniverse(2, 2, fill)
    call map % add(2, 2)

    fill = [3, 2, 1]
    call graph % addUniverse(3, 3, fill)
    call map % add(3, 3)

    call graph % finishBuild(map)
    call graph % setRoot(1)

    @assertTrue(graph % hasCycles())

    ! Move recursion to 3rd diffrent universe (loop from 3rd to 2nd)
    !
    graph % uni(2) % fill(1) = 3
    graph % uni(3) % fill(1) = -2
    @assertTrue(graph % hasCycles())

    ! Remove recursion
    graph % uni(3) % fill(1) = 3

    @assertFalse(graph % hasCycles())

    ! Clean up
    call graph % kill()

  end subroutine test_cycles

  !!
  !! Test counting nesting depth
  !!
@Test
  subroutine test_nesting_count()

    ! Verify depth
    @assertEqual(3, geom % maxNesting())

    ! Decrease to single level
    geom % uni(1) % fill(1) = 8
    @assertEqual(1, geom % maxNesting())

  end subroutine test_nesting_count

  !!
  !! Test search for outside
  !!
@Test
  subroutine test_outside_search()

    @assertTrue(geom % nestedOutside())

    ! Remove the outside
    geom % uni(3) % fill(2) = 8
    @assertFalse(geom % nestedOutside())

  end subroutine test_outside_search

  !!
  !! Find unused universes
  !!
@Test
  subroutine test_unused_universes()
    integer(shortInt), dimension(:), allocatable :: unused
    integer(shortInt)                            :: pos

    unused = geom % unusedUniverses()

    ! Search
    pos = linFind(unused, 200)
    @assertTrue(pos /= targetNotFound, 'Missing ID in unused universes')

    pos = linFind(unused, 201)
    @assertTrue(pos /= targetNotFound, 'Missing ID in unused universes')

  end subroutine test_unused_universes

  !!
  !! Test instances count
  !!
@Test
  subroutine test_count_instances()
    type(intMap) :: map
    integer(shortInt) :: idx

    ! Perform count
    call geom % countInstances(map)

    ! Verify absent -> ID 200 & 2001
    @assertEqual(0, map % get(6))
    @assertEqual(0, map % get(7))

    ! Verfiy multiple instances -> ID 1001
    @assertEqual(2, map % get(3))

    ! Verify single instance
    @assertEqual(1, map % get(1))
    @assertEqual(1, map % get(2))

  end subroutine test_count_instances

end module uniFills_test
