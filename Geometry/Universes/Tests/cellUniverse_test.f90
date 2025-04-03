module cellUniverse_test

  use numPrecision
  use genericProcedures
  use universalVariables, only : UNDEF_MAT
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use charMap_class,      only : charMap
  use coord_class,        only : coord
  use surfaceShelf_class, only : surfaceShelf
  use cellShelf_class,    only : cellShelf
  use cellUniverse_class, only : cellUniverse
  use funit

  implicit none
  
  ! Parameters
  character(*), parameter :: SURF_DEF = &
  " surf1 { id 1; type sphere; origin (0.0 0.0 0.0); radius 2;}&
  & surf2 { id 2; type sphere; origin (4.0 0.0 0.0); radius 1;}"

  character(*), parameter :: CELL_DEF = &
  " cell1 {id 1; type simpleCell; surfaces (-1); filltype uni; universe 3;} &
  & cell2 {id 2; type simpleCell; surfaces (1 2); filltype uni; universe 4;}"


  !
  ! Note that rotation is such that following axis transformation applies:
  !   x -> z
  !   y -> -y
  !   z -> x
  !
  character(*), parameter :: UNI_DEF = &
  "id 1; type cellUniverse; origin (2.0 0.0 0.0); rotation (90.0 90.0 90.0); cells (1 2);"

  ! Variables
  type(surfaceShelf) :: surfs
  type(cellShelf)    :: cells
  type(charMap)      :: mats
  type(cellUniverse) :: uni

contains

  !!
  !! Setup environment
  !!
@Before
  subroutine setUp()
    integer(shortInt), dimension(:), allocatable :: fill
    type(dictionary) :: dict

    ! Build surfaces and MATS

    call charToDict(dict, SURF_DEF)
    call surfs % init(dict)
    call dict % kill()

    call charToDict(dict, CELL_DEF)
    call cells % init(dict, surfs, mats)
    call dict % kill()

    ! Build universe
    call charToDict(dict, UNI_DEF)
    call uni % init(fill, dict, cells, surfs, mats)
    call dict % kill()

    ! Set index
    call uni % setIdx(8)

    ! Verify fill
    @assertEqual([-3, -4, UNDEF_MAT], fill)

  end subroutine setUp

  !!
  !! Clean environment
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

    ! Get id
    @assertEqual(1, uni % id())

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
    r = [0.0_defReal, 0.0_defReal, 3.0_defReal ]
    dir = [ZERO, ZERO, ONE]

    call uni % enter(new, r, dir)

    ! Verify location
    r_ref = [1.0_defReal, 0.0_defReal, 0.0_defReal]
    u_ref = [ONE, ZERO, ZERO]
    @assertEqual(r_ref, new % r, TOL )
    @assertEqual(u_ref, new % dir, TOL)
    @assertEqual(8, new % uniIdx)
    @assertEqual(1, new % localID)
    @assertEqual(cells % getIdx(1), new % cellIdx)

    ! ** Enter into local cell 2
    r = [2.0_defReal, 0.0_defReal, 1.0_defReal]
    dir = [ZERO, ONE, ZERO]

    call uni % enter(new, r, dir)

    ! Verify location
    r_ref = [-1.0_defReal, 0.0_defReal, 2.0_defReal]
    u_ref = [ZERO, -ONE, ZERO]
    @assertEqual(r_ref, new % r, TOL )
    @assertEqual(u_ref, new % dir, TOL)
    @assertEqual(8, new % uniIdx)
    @assertEqual(2, new % localID)
    @assertEqual(cells % getIdx(2), new % cellIdx)

    ! ** Enter into the UNDEFINED cell
    r = [0.0_defReal, 0.0_defReal, 6.5_defReal]
    dir = [ONE, ZERO, ZERO]

    call uni % enter(new, r, dir)

    ! Verify location
    r_ref = [4.5_defReal, 0.0_defReal, 0.0_defReal]
    u_ref = [ZERO, ZERO, ONE]
    @assertEqual(r_ref, new % r, TOL )
    @assertEqual(u_ref, new % dir, TOL)
    @assertEqual(8, new % uniIdx)
    @assertEqual(3, new % localID)
    @assertEqual(0, new % cellIdx)

    ! Verify rotation settings in coord
    ! * Do it only once
    @assertTrue(new % isRotated)
    @assertEqual([ZERO, ZERO,  ONE], new % rotMat(1,:), TOL)
    @assertEqual([ZERO, -ONE, ZERO], new % rotMat(2,:), TOL)
    @assertEqual([ONE , ZERO, ZERO], new % rotMat(3,:), TOL)

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
    pos % r = [-1.0_defReal, 0.0_defReal, 0.0_defReal]
    pos % dir = [ONE, ZERO, ZERO]
    pos % uniIdx  = 8
    pos % cellIdx = cells % getIdx(1)
    pos % localId = 1

    call uni % distance(d, surfIdx, pos)

    ref = 3.0_defReal
    @assertEqual(ref, d, TOL * ref)
    @assertEqual(surfs % getIdx(1), surfIdx )


    ! ** In local cell 2 distance to surface 2
    pos % r = [7.0_defReal, 0.0_defReal, 0.0_defReal]
    pos % dir = [-ONE, ZERO, ZERO]
    pos % cellIdx = cells % getIdx(2)
    pos % localId = 2

    call uni % distance(d, surfIdx, pos)

    ref = 2.0_defReal
    @assertEqual(ref, d, TOL * ref)
    @assertEqual(surfs % getIdx(2), surfIdx )

    ! ** In local cell 2 distance to infinity
    ! surfIdx must be set to 0
    pos % dir = [ONE, ZERO, ZERO]
    call uni % distance(d, surfIdx, pos)

    @assertEqual(INF, d)
    @assertEqual(0, surfIdx)

  end subroutine test_distance

  !!
  !! Test cell-to cell crossing
  !!
@Test
  subroutine test_cross()
    type(coord)       :: pos
    integer(shortInt) :: idx

    ! Cross from cell 1 to cell 2
    pos % r   = [0.0_defReal, 2.0_defReal, 0.0_defReal]
    pos % dir = [ZERO, ONE, ZERO]
    pos % uniIdx = 8
    pos % cellIdx = cells % getIdx(1)
    pos % localId = 1

    idx = surfs % getIdx(1)
    call uni % cross(pos, idx)

    @assertEqual(2, pos % localId)
    @assertEqual(cells % getIdx(2), pos % cellIdx)

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
    pos % uniIdx = 8
    pos % cellIdx = cells % getIdx(1)
    pos % localId = 1

    @assertEqual([ZERO, ZERO, ZERO], uni % cellOffset(pos) )

    ! Cell 2
    pos % r   = [-7.0_defReal, 2.0_defReal, 0.0_defReal]
    pos % dir = [ZERO, ONE, ZERO]
    pos % uniIdx = 8
    pos % cellIdx = cells % getIdx(2)
    pos % localId = 2

    @assertEqual([ZERO, ZERO, ZERO], uni % cellOffset(pos) )

  end subroutine test_cellOffset



end module cellUniverse_test
