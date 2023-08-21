module pinUniverse_test

  use numPrecision
  use universalVariables, only : INF, SURF_TOL
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use charMap_class,      only : charMap
  use coord_class,        only : coord
  use surfaceShelf_class, only : surfaceShelf
  use cellShelf_class,    only : cellShelf
  use pinUniverse_class,  only : pinUniverse, MOVING_IN, MOVING_OUT
  use funit
  implicit none

  ! Parameters
  character(*), parameter :: UNI_DEF = &
  "id 7; type pinUniverse; origin (0.0 0.0 0.0); rotation (0.0 0.0 0.0); &
  &radii (2.5 1.5 0.0); fills (u<7> u<14> void);"

  ! Variables
  type(surfaceShelf) :: surfs
  type(cellShelf)    :: cells
  type(charMap)      :: mats
  type(pinUniverse)  :: uni


contains

  !!
  !! Set-up test enviroment
  !!
@Before
  subroutine setup()
    character(nameLen)                           :: name
    integer(shortInt), dimension(:), allocatable :: fill
    type(dictionary)                             :: dict

    ! Load void material
    name = 'void'
    call mats % add(name, 13)

    ! Build universe
    call charToDict(dict, UNI_DEF)
    call uni % init(fill, dict, cells, surfs, mats)

    ! Set index
    call uni % setIdx(3)

    ! Verify fill array
    @assertEqual([-14, -7, 13], fill)


  end subroutine setup

  !!
  !! Clean after test
  !!
@After
  subroutine clean()

    call surfs % kill()
    call cells % kill()
    call mats % kill()
    call uni % kill()

  end subroutine clean

  !!
  !! Test miscellaneous functionality
  !!
@Test
  subroutine test_misc()

    ! Get id
    @assertEqual(7, uni % id())

    ! Set ID
    call uni % setId(7)
    @assertEqual(7, uni % id())

  end subroutine test_misc

  !!
  !! Test entering a universe
  !!
@Test
  subroutine test_enter()
    type(coord) :: new
    real(defReal), dimension(3) :: r_ref, u_ref, r, dir
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    ! ** Enter into local cell 1
    r = [0.0_defReal, 1.0_defReal, 0.0_defReal ]
    dir = [ZERO, ZERO, ONE]

    call uni % enter(new, r, dir)

    ! Verify location
    r_ref = r
    u_ref = dir
    @assertEqual(r_ref, new % r, TOL)
    @assertEqual(u_ref, new % dir, TOL)
    @assertEqual(3, new % uniIdx)
    @assertEqual(1, new % localID)
    @assertEqual(0, new % cellIdx)

    ! ** Enter into local cell 2
    r = [2.3_defReal, 0.0_defReal, -980.0_defReal ]
    dir = [ZERO, ZERO, ONE]

    call uni % enter(new, r, dir)

    ! Verify location
    r_ref = r
    u_ref = dir
    @assertEqual(r_ref, new % r, TOL)
    @assertEqual(u_ref, new % dir, TOL)
    @assertEqual(3, new % uniIdx)
    @assertEqual(2, new % localID)
    @assertEqual(0, new % cellIdx)

    ! ** Enter into local cell 3
    r = [2.6_defReal, 0.0_defReal, -980.0_defReal ]
    dir = [ZERO, ZERO, ONE]

    call uni % enter(new, r, dir)

    ! Verify location
    r_ref = r
    u_ref = dir
    @assertEqual(r_ref, new % r, TOL)
    @assertEqual(u_ref, new % dir, TOL)
    @assertEqual(3, new % uniIdx)
    @assertEqual(3, new % localID)
    @assertEqual(0, new % cellIdx)

    ! VERIFY THAT ROTATION IS NOT SET (all angles were 0.0)
    @assertFalse(new % isRotated)

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

    ! ** In local cell 1 distance to boundary
    pos % r = [1.0_defReal, 0.0_defReal, 0.0_defReal]
    pos % dir = [ONE, ZERO, ZERO]
    pos % uniIdx  = 3
    pos % cellIdx = 0
    pos % localId = 1

    call uni % distance(d, surfIdx, pos)

    ref = 0.5_defReal
    @assertEqual(ref, d, ref * tol)
    @assertEqual(MOVING_OUT, surfIdx)

    ! ** In outermost cell moving away
    pos % r = [2.0_defReal, 1.6_defReal, 0.0_defReal]
    pos % dir = [ONE, ZERO, ZERO]
    pos % localId = 3

    call uni % distance(d, surfIdx, pos)
    @assertEqual(INF, d)
    ! Surface momento is undefined -> No crossing

    ! In ordinary cell in-between
    pos % r = [0.0_defReal, 1.6_defReal, 0.0_defReal]
    pos % dir = [ZERO, -ONE, ZERO]
    pos % localId = 2

    call uni % distance(d, surfIdx, pos)
    ref = 0.1_defReal
    @assertEqual(ref, d, ref * tol)
    @assertEqual(MOVING_IN, surfIdx)

  end subroutine test_distance

  !!
  !! Test cell-to cell crossing
  !!
@Test
  subroutine test_cross()
    type(coord)       :: pos
    integer(shortInt) :: idx
    real(defReal) :: eps

    ! Cross from cell 1 to cell 2
    eps = HALF * SURF_TOL
    pos % r   = [0.0_defReal, 1.5_defReal-eps, 0.0_defReal]
    pos % dir = [ZERO, ONE, ZERO]
    pos % uniIdx = 8
    pos % cellIdx = 0
    pos % localId = 1

    idx = MOVING_OUT
    call uni % cross(pos, idx)

    @assertEqual(2, pos % localId)

    ! Cross form cell 2 to cell 1
    eps = HALF * SURF_TOL
    pos % r   = [0.0_defReal, 1.5_defReal+eps, 0.0_defReal]
    pos % dir = [ZERO, -ONE, ZERO]

    idx = MOVING_IN
    call uni % cross(pos, idx)

    @assertEqual(1, pos % localId)

  end subroutine test_cross

  !!
  !! Test cell offset
  !!
@Test
  subroutine test_cellOffset()
    type(coord)       :: pos

    ! Cell 1
    pos % r   = [0.0_defReal, 1.0_defReal, 0.0_defReal]
    pos % dir = [ZERO, ONE, ZERO]
    pos % uniIdx = 3
    pos % cellIdx = 0
    pos % localId = 1

    @assertEqual([ZERO, ZERO, ZERO], uni % cellOffset(pos) )

    ! Cell 3
    pos % r   = [-7.0_defReal, 2.0_defReal, 0.0_defReal]
    pos % dir = [ZERO, ONE, ZERO]
    pos % uniIdx = 3
    pos % cellIdx = 0
    pos % localId = 3

    @assertEqual([ZERO, ZERO, ZERO], uni % cellOffset(pos) )

  end subroutine test_cellOffset

  !!
  !! Test surface transitions
  !!
  !! Check that there is no problem with distance calculations
  !! if particle is placed very close to an annulus surface (within SURF_TOL)
  !!
@Test
  subroutine test_edgeCases()
    type(coord)       :: pos
    integer(shortInt) :: idx, localID, cellIdx
    real(defReal)     :: eps, d
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    ! At boundary between cell 1 and 2
    eps = HALF * SURF_TOL
    pos % r   = [0.0_defReal, 1.5_defReal-eps, 0.0_defReal]
    pos % dir = [ONE, -0.00001_defReal, ZERO]
    pos % dir = pos % dir / norm2(pos % dir)
    pos % uniIdx = 8
    pos % cellIdx = 0

    ! Should find particle in cell 1
    ! And return very small distance -> MOVING OUT
    call uni % findCell(localID, cellIDx, pos % r, pos % dir)
    @assertEqual(1, localID)

    pos % localID = 1
    call uni % distance(d, idx, pos)

    @assertEqual(ZERO, d, 1.0E-3_defReal)
    @assertEqual(MOVING_OUT, idx)

  end subroutine test_edgeCases


end module pinUniverse_test
