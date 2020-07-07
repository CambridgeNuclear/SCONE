module rootUniverse_test

  use numPrecision
  use universalVariables, only : OUTSIDE_MAT
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use charMap_class,      only : charMap
  use coord_class,        only : coord
  use surfaceShelf_class, only : surfaceShelf
  use cellShelf_class,    only : cellShelf
  use rootUniverse_class, only : rootUniverse
  use pfUnit_mod

  implicit none

  ! Parameters
  character(*), parameter :: SURF_DEF = &
  "surf1 { id 4; type sphere; origin (0.0 5.0 0.0); radius 0.5;}&
  &surf2 { id 1; type sphere; origin (0.0 0.0 0.0); radius 2;}"

  character(*), parameter :: UNI_DEF = &
  "id 1; type rootUniverse; border 1; fill u<17>;"

  ! Variables
  type(surfaceShelf) :: surfs
  type(cellShelf)    :: cells
  type(charMap)      :: mats
  type(rootUniverse) :: uni

contains

  !!
  !! Setup enviroment
  !!
@Before
  subroutine setUp()
    integer(shortInt), dimension(:), allocatable :: fill
    type(dictionary) :: dict

    ! Build surfaces and MATS

    call charToDict(dict, SURF_DEF)
    call surfs % init(dict)
    call dict % kill()

    ! Build universe
    call charToDict(dict, UNI_DEF)
    call uni % init(fill, dict, cells, surfs, mats)
    call dict % kill()

    ! Set index
    call uni % setIdx(8)

    ! Verify fill
    @assertEqual([-17, OUTSIDE_MAT], fill)

  end subroutine setUp

  !!
  !! Clean enviroment
  !!
@After
  subroutine clean()

    call surfs % kill()
    call cells % kill()
    call mats % kill()
    call uni % kill()

  end subroutine clean

  !!
  !! Test miscellaneous functionality (of generic universe)
  !!
@Test
  subroutine test_misc()
    real(defReal), dimension(3,3) :: mat

    ! Get id
    @assertEqual(1, uni % id())

    ! Set ID
    call uni % setId(7)
    @assertEqual(7, uni % id())

    ! Test boundary surface
    @assertEqual(surfs % getIdx(1), uni % border())

  end subroutine test_misc

  !!
  !! Test entering a universe
  !!
@Test
  subroutine test_enter()
    type(coord) :: old
    type(coord) :: new
    real(defReal), dimension(3) :: r_ref, u_ref
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    ! Enter inside
    old % r = [1.0_defReal, -1.0_defReal, 1.0_defReal]
    old % dir = [ONE, ZERO, ZERO]

    call uni % enter(new, old)

    @assertEqual(old % r, new % r, TOL)
    @assertEqual(old % dir, new % dir, TOL)
    @assertEqual(8, new % uniIdx)
    @assertEqual(0, new % cellIdx)
    @assertEqual(1, new % localID)

    ! Enter outside
    old % r = [2.0_defReal, -2.0_defReal, 1.0_defReal]
    old % dir = [ONE, ZERO, ZERO]

    call uni % enter(new, old)

    @assertEqual(old % r, new % r, TOL)
    @assertEqual(old % dir, new % dir, TOL)
    @assertEqual(8, new % uniIdx)
    @assertEqual(0, new % cellIdx)
    @assertEqual(2, new % localID)

  end subroutine test_enter

  !!
  !! Test distance calculation
  !!
@Test
  subroutine test_distance()
    real(defReal)     :: d, ref
    integer(shortInt) :: surfIdx
    type(coord)       :: pos
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    ! Distance from inside -> only relevant
    pos % r = [1.0_defReal, 0.0_defReal, 0.0_defReal]
    pos % dir = [ONE, ZERO, ZERO]
    pos % uniIdx  = 8
    pos % localID = 1

    call uni % distance(d, surfIdx, pos)

    ref = 1.0_defReal
    @assertEqual(ref, d, ref * TOL)
    @assertEqual(surfs % getIdx(1), surfIdx)

  end subroutine test_distance

  !!
  !! Test cell-to cell crossing
  !!
@Test
  subroutine test_cross()
    type(coord)       :: pos
    integer(shortInt) :: idx

    ! Cross into outside
    pos % r = [2.0_defReal, 0.0_defReal, 0.0_defReal]
    pos % dir = [ONE, ZERO, ZERO]
    pos % uniIdx  = 8
    pos % localID = 1

    idx = surfs % getIdx(1)
    call uni % cross(pos, idx)

    @assertEqual(2, pos % localID)

  end subroutine test_cross

  !!
  !! Test cell offset
  !!
@Test
  subroutine test_cellOffset()
    type(coord)       :: pos

    ! Inside
    pos % r = [1.5_defReal, 0.0_defReal, 0.0_defReal]
    pos % dir = [ONE, ZERO, ZERO]
    pos % uniIdx  = 8
    pos % localID = 1

    @assertEqual([ZERO, ZERO, ZERO], uni % cellOffset(pos))

    ! Outside
    pos % r = [2.5_defReal, 0.0_defReal, 0.0_defReal]
    pos % dir = [ONE, ZERO, ZERO]
    pos % uniIdx  = 8
    pos % localID = 2

    @assertEqual([ZERO, ZERO, ZERO], uni % cellOffset(pos))

  end subroutine test_cellOffset



end module rootUniverse_test
