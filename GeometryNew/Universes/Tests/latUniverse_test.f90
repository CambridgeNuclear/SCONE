module latUniverse_test

  use numPrecision
  use genericProcedures
  use universalVariables, only : UNDEF_MAT
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use charMap_class,      only : charMap
  use coord_class,        only : coord
  use surfaceShelf_class, only : surfaceShelf
  use cellShelf_class,    only : cellShelf
  use latUniverse_class,  only : latUniverse
  use pfUnit_mod

  implicit none

  ! Parameters
  character(*), parameter :: UNI1_DEF = &
  "id 1; type latUniverse; origin (0.0 0.0 0.0); rotation (0.0 0.0 0.0); &
   pitch (1.0 2.0 3.0); shape (3 2 2); padMat void; &
  &map ( 3 4 5 &
  &      7 4 8 &
  &            &
  &      1 2 3 &
  &      4 5 6); "

  character(*), parameter :: UNI2_DEF = &
  "id 2; type latUniverse; pitch (1.0 2.0 0.0); shape (2 1 0); padMat u<1>; &
  &map (1 2); "

  ! Variables
  type(surfaceShelf) :: surfs
  type(cellShelf)    :: cells
  type(charMap)      :: mats
  type(latUniverse)  :: uni1
  type(latUniverse)  :: uni2

contains

  !!
  !! Setup enviroment
  !!
@Before
  subroutine setUp()
    integer(shortInt), dimension(:), allocatable :: fill
    type(dictionary)   :: dict
    character(nameLen) :: name
    integer(shortInt), dimension(:), allocatable :: ref

    ! Add materials
    name = 'void'
    call mats % add(name, 3)

    ! Build universe 1
    call charToDict(dict, UNI1_DEF)
    call uni1 % init(fill, dict, cells, surfs, mats)
    call dict % kill()
    call uni1 % setIdx(8)

    ! Verify fill vector
    ref = [-4, -5, -6, -1, -2, -3, -7, -4, -8, -3, -4, -5, 3]
    @assertEqual(ref, fill)

    ! Build universe 2
    call charToDict(dict, UNI2_DEF)
    call uni2 % init(fill, dict, cells, surfs, mats)
    call dict % kill()
    call uni2 % setIdx(3)

    ref = [-1, -2, -1]
    @assertEqual(ref, fill)

  end subroutine setUp

  !!
  !! Clean enviroment
  !!
@After
  subroutine clean()

    call surfs % kill()
    call cells % kill()
    call mats % kill()
    call uni1 % kill()
    call uni2 % kill()

  end subroutine clean

  !!
  !! Test miscellaneous functionality (of generic universe)
  !!
@Test
  subroutine test_misc()
    real(defReal), dimension(3,3) :: mat

    ! * Single universe is fine here
    ! Get id
    @assertEqual(1, uni1 % id())

    ! Set ID
    call uni1 % setId(7)
    @assertEqual(7, uni1 % id())

  end subroutine test_misc

  !!
  !! Test entering a universe
  !!
@Test
  subroutine test_enter()
    type(coord) :: new
    real(defReal), dimension(3) :: r_ref, u_ref, r ,dir
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    ! ** 3D universe
    ! Enter inside -> Away from surface
    r = [1.0_defReal, 1.0_defReal, 0.5_defReal ]
    dir = [ZERO, ZERO, ONE]

    call uni1 % enter(new, r, dir)

    r_ref = r
    u_ref = dir
    @assertEqual(r_ref, new % r, TOL )
    @assertEqual(u_ref, new % dir, TOL)
    @assertEqual(8, new % uniIdx)
    @assertEqual(12, new % localID)
    @assertEqual(0, new % cellIdx)

    ! Enter outside
    r = [1.6_defReal, 0.5_defReal, 0.5_defReal ]
    dir = [ZERO, ZERO, ONE]

    call uni1 % enter(new, r, dir)

    r_ref = r
    u_ref = dir
    @assertEqual(r_ref, new % r, TOL )
    @assertEqual(u_ref, new % dir, TOL)
    @assertEqual(8, new % uniIdx)
    @assertEqual(13, new % localID)
    @assertEqual(0, new % cellIdx)

    ! Enter in a corner
    r = [-0.5_defReal, 0.0_defReal, 0.0_defReal ]
    dir = [-ONE, ONE, -ONE]
    dir = dir / norm2(dir)

    call uni1 % enter(new, r, dir)

    r_ref = r
    u_ref = dir
    @assertEqual(r_ref, new % r, TOL )
    @assertEqual(u_ref, new % dir, TOL)
    @assertEqual(8, new % uniIdx)
    @assertEqual(4, new % localID)
    @assertEqual(0, new % cellIdx)

    ! ** 2D Universe
    ! Enter inside -> Away from surface
    r = [0.5_defReal, 0.5_defReal, 13.5_defReal ]
    dir = [ZERO, ZERO, ONE]

    call uni2 % enter(new, r, dir)

    r_ref = r
    u_ref = dir
    @assertEqual(r_ref, new % r, TOL )
    @assertEqual(u_ref, new % dir, TOL)
    @assertEqual(3, new % uniIdx)
    @assertEqual(2, new % localID)
    @assertEqual(0, new % cellIdx)

    ! Enter outside
    r = [1.6_defReal, 0.5_defReal, 0.5_defReal ]
    dir = [ZERO, ZERO, ONE]

    call uni2 % enter(new, r, dir)

    r_ref = r
    u_ref = dir
    @assertEqual(r_ref, new % r, TOL )
    @assertEqual(u_ref, new % dir, TOL)
    @assertEqual(3, new % uniIdx)
    @assertEqual(3, new % localID)
    @assertEqual(0, new % cellIdx)

    ! Enter on a face
    r = [0.0_defReal, 0.0_defReal, 0.0_defReal ]
    dir = [-ONE, ONE, -ONE]
    dir = dir / norm2(dir)

    call uni2 % enter(new, r, dir)

    r_ref = r
    u_ref = dir
    @assertEqual(r_ref, new % r, TOL )
    @assertEqual(u_ref, new % dir, TOL)
    @assertEqual(3, new % uniIdx)
    @assertEqual(1, new % localID)
    @assertEqual(0, new % cellIdx)


  end subroutine test_enter

  !!
  !! Test distance calculation
  !!
@Test
  subroutine test_distance()
    real(defReal)     :: d, ref, eps
    integer(shortInt) :: surfIdx
    type(coord)       :: pos
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    ! ** 3D universe
    ! Well inside a cell
    pos % r = [0.0_defReal, 0.1_defReal, 0.5_defReal]
    pos % dir = [ZERO, -ZERO, ONE]
    pos % uniIdx  = 8
    pos % cellIdx = 0
    pos % localId = 11

    call uni1 % distance(d, surfIdx, pos)

    ref = 2.5_defReal
    @assertEqual(ref, d, TOL * ref)
    @assertEqual(-6, surfIdx )

    ! From outside -> miss
    pos % r = [-4.0_defReal, 0.1_defReal, 0.5_defReal]
    pos % dir = [ONE, ONE, ZERO]
    pos % dir = pos % dir / norm2(pos % dir)
    pos % uniIdx  = 8
    pos % cellIdx = 0
    pos % localId = 13

    call uni1 % distance(d, surfIdx, pos)

    @assertEqual(INF, d)
    @assertEqual(-7, surfIdx )

    ! After a surface undershoot
    eps = HALF * SURF_TOL
    pos % r = [-1.0_defReal, 0.0_defReal-eps, -0.5_defReal]
    pos % dir = [ONE, ONE, ZERO]
    pos % dir = pos % dir / norm2(pos % dir)
    pos % uniIdx  = 8
    pos % cellIdx = 0
    pos % localId = 4

    call uni1 % distance(d, surfIdx, pos)

    ref = SQRT2 * HALF
    @assertEqual(ref, d, ref * TOL)
    @assertEqual(-2, surfIdx )

    ! After overshoot via a corner
    pos % r = [-0.5_defReal+eps, 0.0_defReal+eps, -0.5_defReal]
    pos % dir = [ONE, ONE, ZERO]
    pos % dir = pos % dir / norm2(pos % dir)
    pos % uniIdx  = 8
    pos % cellIdx = 0
    pos % localId = 4

    call uni1 % distance(d, surfIdx, pos)
    @assertEqual(ZERO, d,  TOL)
    @assertEqual(-2, surfIdx )

    !** 2D universe
    ! Well inside a cell -> Vertical
    pos % r = [0.5_defReal, 0.6_defReal, 0.5_defReal]
    pos % dir = [ZERO, ZERO, ONE]
    pos % uniIdx  = 3
    pos % cellIdx = 0
    pos % localId = 2

    call uni2 % distance(d, surfIdx, pos)

    @assertEqual(INF, d)

    ! Well inside a cell -> Shallow hit
    pos % r = [0.5_defReal, 0.6_defReal, 0.5_defReal]
    pos % dir = [ZERO, 0.01_defReal, ONE]
    pos % dir = pos % dir / norm2(pos % dir)

    call uni2 % distance(d, surfIdx, pos)

    ref = sqrt(40.0_defReal**2 + 0.4_defReal**2)
    @assertEqual(ref, d, TOL * ref)
    @assertEqual(-4, surfIdx )

    ! From outside -> Hit
    pos % r = [-1.5_defReal, 0.6_defReal, 0.5_defReal]
    pos % dir = [ONE, ZERO, ONE]
    pos % dir = pos % dir / norm2(pos % dir)
    pos % uniIdx  = 3
    pos % cellIdx = 0
    pos % localId = 3

    call uni2 % distance(d, surfIdx, pos)

    ref = HALF * SQRT2
    @assertEqual(ref, d, TOL * ref)
    @assertEqual(-7, surfIdx )

  end subroutine test_distance

  !!
  !! Test cell-to cell crossing
  !!
@Test
  subroutine test_cross()
    type(coord)       :: pos
    integer(shortInt) :: idx

    ! *** 3D Lattice
    ! Cross inside
    pos % r = [-1.0_defReal, 0.0_defReal, 0.5_defReal]
    pos % dir = [-ONE, ONE, -ONE]
    pos % dir = pos % dir / norm2(pos % dir)
    pos % uniIdx  = 8
    pos % cellIdx = 0
    pos % localId = 1

    call uni1 % cross(pos, -4)

    @assertEqual(4, pos % localID)

    ! Cross from outside
    pos % r = [1.0_defReal, 2.0_defReal, -0.5_defReal]
    pos % dir = [ONE, -ONE, -ONE]
    pos % dir = pos % dir / norm2(pos % dir)
    pos % localId = 13

    call uni1 % cross(pos, -7)

    @assertEqual(6, pos % localID)

    ! *** 2D Lattice
    pos % r = [0.0_defReal, 0.0_defReal, 16.5_defReal]
    pos % dir = [ONE, ONE, -ONE]
    pos % dir = pos % dir / norm2(pos % dir)
    pos % uniIdx  = 3
    pos % cellIdx = 0
    pos % localId = 1

    call uni2 % cross(pos, -2)

    @assertEqual(2, pos % localID)

    ! Cross from outside
    pos % r = [-1.0_defReal, -0.5_defReal, -78.5_defReal]
    pos % dir = [ONE, ONE, ZERO]
    pos % dir = pos % dir / norm2(pos % dir)
    pos % localId = 3

    call uni2 % cross(pos, -7)

    @assertEqual(1, pos % localID)

  end subroutine test_cross

  !!
  !! Test cell offset
  !!
@Test
  subroutine test_cellOffset()
    type(coord)                 :: pos
    real(defReal), dimension(3) :: ref
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    ! ** 3D lattice
    ! Inside
    pos % r = [0.0_defReal, 0.0_defReal, 0.5_defReal]
    pos % dir = [-ONE, ONE, -ONE]
    pos % dir = pos % dir / norm2(pos % dir)
    pos % uniIdx  = 8
    pos % cellIdx = 0
    pos % localId = 11

    ref = [0.0_defReal, 1.0_defReal, 1.5_defReal]
    @assertEqual(ref, uni1 % cellOffset(pos), TOL)

    ! Outside
    pos % r = [-7.0_defReal, 0.0_defReal, 0.5_defReal]
    pos % dir = [-ONE, ONE, -ONE]
    pos % dir = pos % dir / norm2(pos % dir)
    pos % localId = 13

    ref = ZERO
    @assertEqual(ref, uni1 % cellOffset(pos), TOL)

    ! ** 2D Lattice
    pos % r = [0.5_defReal, 0.0_defReal, 0.5_defReal]
    pos % dir = [-ONE, ONE, -ONE]
    pos % dir = pos % dir / norm2(pos % dir)
    pos % uniIdx  = 3
    pos % cellIdx = 0
    pos % localId = 2

    ref = [0.5_defReal, 0.0_defReal, 0.0_defReal]
    @assertEqual(ref, uni2 % cellOffset(pos), TOL)

    ! Outside
    pos % r = [-7.0_defReal, 0.0_defReal, 0.5_defReal]
    pos % dir = [-ONE, ONE, -ONE]
    pos % dir = pos % dir / norm2(pos % dir)
    pos % localId = 3

    ref = ZERO
    @assertEqual(ref, uni2 % cellOffset(pos), TOL)

  end subroutine test_cellOffset


end module latUniverse_test
