module cellShelf_test
  use numPrecision
  use universalVariables, only : OUTSIDE_MAT, VOID_MAT
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use charMap_class,      only : charMap
  use surfaceShelf_class, only : surfaceShelf
  use cellShelf_class,    only : cellShelf
  use cell_inter,         only : cell
  use funit

  implicit none

  ! Parameters
  character(*), parameter :: SURFS_DEF = "&
  & surf1 {id 7; type sphere; origin (0.0 0.0 0.0); radius 2.0;} &
  & surf2 {id 8; type zPlane; z0 -1.3;} &
  & surf3 {id 9; type zPlane; z0 0.0;} &
  & surf4 {id 10; type zPlane; z0 1.0;}"

  character(*), parameter :: CELLS_DEF = &
  " cell1 {id 4; type simpleCell; surfaces (7); filltype outside;} &
  & cell2 {id 2; type simpleCell; surfaces (-7 10); filltype uni; universe 8;} &
  & cell3 {id 1; type simpleCell; surfaces (-7 -10 9); filltype mat; material void;} &
  & cell4 {id 5; type simpleCell; surfaces (-7 -9 8); filltype mat; material fuel;} &
  & cell5 {id 8; type simpleCell; surfaces (-7 -8); filltype mat; material pecorino;}"


  ! Variables
  type(charMap)      :: mats
  type(surfaceShelf) :: surfs
  type(cellShelf)    :: cells

contains

  !!
  !! Build the cell
  !!
@Before
  subroutine setUp()
    type(dictionary) :: dict
    character(nameLen) :: name

    ! Surface shelf
    call charToDict(dict, SURFS_DEF)
    call surfs % init(dict)
    call dict % kill()

    ! Material map
    name = 'void'
    call mats % add(name, VOID_MAT)
    name = 'fuel'
    call mats % add(name, 1)
    name = 'pecorino'
    call mats % add(name, 2)

    ! Build cells
    call charToDict(dict, CELLS_DEF)
    call cells % init(dict, surfs, mats)

  end subroutine setUp

  !!
  !! Clean after tests
  !!
@After
  subroutine cleanUp()

    call surfs % kill()
    call cells % kill()
    call mats % kill()

  end subroutine cleanUp

  !!
  !! Test shelf
  !!
@Test
  subroutine test_get()
    integer(shortInt) :: idx
    class(cell), pointer :: ptr

    ! Cell ID 1
    idx = cells % getIdx(1)
    ptr => cells % getPtr(idx)
    @assertEqual(1, ptr % id())
    @assertEqual(1, cells % getID(idx))

    ! Cell ID 5
    idx = cells % getIdx(5)
    ptr => cells % getPtr(idx)
    @assertEqual(5, ptr % id())
    @assertEqual(5, cells % getID(idx))

    ! Cell ID 8
    idx = cells % getIdx(8)
    ptr => cells % getPtr(idx)
    @assertEqual(8, ptr % id())
    @assertEqual(8, cells % getID(idx))

    ! Test size
    @assertEqual(5, cells % getSize())

  end subroutine test_get

  !!
  !! Test cell filling
  !!
@Test
  subroutine test_cell_fill()
    integer(shortInt) :: idx

    ! Cell ID 4
    idx = cells % getIdx(4)
    @assertEqual(OUTSIDE_MAT, cells % getFill(idx))

    ! Cell ID 2 -> universe
    idx = cells % getIdx(2)
    @assertEqual(-8, cells % getFill(idx))

    ! Cell ID 1 -> void
    idx = cells % getIdx(1)
    @assertEqual(VOID_MAT, cells % getFill(idx))

    ! Cell ID 8 -> Pecorino fill
    idx = cells % getIdx(8)
    @assertEqual(2, cells % getFill(idx))

  end subroutine test_cell_fill

end module cellShelf_test
