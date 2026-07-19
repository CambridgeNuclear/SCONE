module hexUniverse_test

  use numPrecision
  use genericProcedures
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use charMap_class,      only : charMap
  use coord_class,        only : coord
  use surfaceShelf_class, only : surfaceShelf
  use cellShelf_class,    only : cellShelf
  use hexUniverse_class
  use funit

  implicit none

  ! Parameters
  ! uni1: 3D lattice, "point" orientation, hexagonal plane in x-y, stacked along z
  !   halfwidth = 1.0 (pitch 2.0), axial pitch = 5.0, shape 3 x 3 x 2
  character(*), parameter :: UNI1_DEF = &
  "id 1; type zHexUniverse; orientation 1; pitch (2.0 5.0); shape (3 3 2); padMat void; &
  & map ( 1 2 3   &
  &       4 5 6   &
  &       7 8 9   &
  &                &
  &      10 11 12 &
  &      13 14 15 &
  &      16 17 18 ); "

  ! uni2: 2D lattice (infinite along the axis), "flat" orientation, hexagonal plane in x-z,
  !   stacked (infinitely) along y. Background filled with a universe.
  character(*), parameter :: UNI2_DEF = &
  "id 2; type yHexUniverse; orientation 2; pitch (2.0 0.0); shape (2 2 0); padMat u<7>; &
  & map ( 1 2   &
  &       2 1); "

  ! Variables
  type(surfaceShelf)   :: surfs
  type(cellShelf)      :: cells
  type(charMap)        :: mats
  type(hexUniverse) :: uni1
  type(hexUniverse) :: uni2

contains

  !!
  !! Setup environment
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
    ref = [-16, -17, -18, -13, -14, -15, -10, -11, -12, &
            -7,  -8,  -9,  -4,  -5,  -6,  -1,  -2,  -3, 3]
    @assertEqual(ref, fill)

    ! Build universe 2
    call charToDict(dict, UNI2_DEF)
    call uni2 % init(fill, dict, cells, surfs, mats)
    call dict % kill()
    call uni2 % setIdx(3)

    ref = [-2, -1, -1, -2, -7]
    @assertEqual(ref, fill)

  end subroutine setUp

  !!
  !! Clean environment
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

    @assertEqual(1, uni1 % id())

    call uni1 % setId(7)
    @assertEqual(7, uni1 % id())

  end subroutine test_misc

  !!
  !! Test entering a universe (findCell through the generic `enter` procedure)
  !!
@Test
  subroutine test_enter()
    type(coord) :: new
    real(defReal), dimension(3) :: r_ref, u_ref, r, dir
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    ! ** 3D lattice
    ! Centre of cell (2,2,1) -> localID 5
    r = [ZERO, ZERO, -2.5_defReal]
    dir = [ZERO, ZERO, ONE]

    call uni1 % enter(new, r, dir)

    r_ref = r
    u_ref = dir
    @assertEqual(r_ref, new % r, TOL)
    @assertEqual(u_ref, new % dir, TOL)
    @assertEqual(8, new % uniIdx)
    @assertEqual(5, new % localID)
    @assertEqual(0, new % cellIdx)

    ! Well outside the lattice -> background
    r = [100.0_defReal, ZERO, -2.5_defReal]
    dir = [ONE, ZERO, ZERO]

    call uni1 % enter(new, r, dir)

    @assertEqual(19, new % localID)

    ! ** 2D (infinite axial) lattice
    ! Centre of cell (2,1) -> localID 2
    r = [0.8660254037844386_defReal, 123.4_defReal, -0.5_defReal]
    dir = [ZERO, ONE, ZERO]

    call uni2 % enter(new, r, dir)

    @assertEqual(3, new % uniIdx)
    @assertEqual(2, new % localID)
    @assertEqual(0, new % cellIdx)

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

    ! ** 3D lattice
    ! Well inside cell 5 (centre), moving along the axis -> hits the axial max face
    pos % r   = [ZERO, ZERO, -2.5_defReal]
    pos % dir = [ZERO, ZERO, ONE]
    pos % uniIdx  = 8
    pos % cellIdx = 0
    pos % localId = 5

    call uni1 % distance(d, surfIdx, pos)

    ref = 2.5_defReal
    @assertEqual(ref, d, TOL * ref)
    @assertEqual(AX_MAX, surfIdx)

    ! Well inside cell 5, moving toward the +x hexagon edge (face 3, negative side)
    pos % r   = [ZERO, ZERO, -2.5_defReal]
    pos % dir = [ONE, ZERO, ZERO]
    pos % localId = 5

    call uni1 % distance(d, surfIdx, pos)

    ref = ONE
    @assertEqual(ref, d, TOL * ref)
    @assertEqual(FACE3_NEG, surfIdx)

    ! Sat exactly on that face, still assigned to cell 5, moving further out -> distance ~ 0
    pos % r   = [ONE, ZERO, -2.5_defReal]
    pos % dir = [ONE, ZERO, ZERO]
    pos % localId = 5

    call uni1 % distance(d, surfIdx, pos)

    @assertEqual(ZERO, d, TOL)
    @assertEqual(FACE3_NEG, surfIdx)

    ! From well outside the lattice, aimed straight at it -> outline surface, then a real hit
    pos % r   = [100.0_defReal, ZERO, -2.5_defReal]
    pos % dir = [-ONE, ZERO, ZERO]
    pos % localId = 19

    call uni1 % distance(d, surfIdx, pos)

    ref = 95.8452995_defReal
    @assertEqual(ref, d, TOL * 10 * ref)
    @assertEqual(OUTLINE_SURF, surfIdx)

  end subroutine test_distance

  !!
  !! Test cell-to-cell crossing
  !!
@Test
  subroutine test_cross()
    type(coord) :: pos

    ! Cross from cell 5 into cell 6 (neighbour across face 3)
    pos % r   = [ONE, ZERO, -2.5_defReal]
    pos % dir = [ONE, ZERO, ZERO]
    pos % uniIdx  = 8
    pos % cellIdx = 0
    pos % localId = 5

    call uni1 % cross(pos, FACE3_NEG)

    @assertEqual(6, pos % localID)

    ! Cross from background back into the lattice (cell 4), at the true tiling boundary
    pos % r   = [-3.0_defReal, ZERO, -2.5_defReal]
    pos % dir = [ONE, ZERO, ZERO]
    pos % localId = 19

    call uni1 % cross(pos, OUTLINE_SURF)

    @assertEqual(4, pos % localID)

    ! Cross out of the lattice, from cell 6 into background
    pos % r   = [3.0_defReal, ZERO, -2.5_defReal]
    pos % dir = [ONE, ZERO, ZERO]
    pos % localId = 6

    call uni1 % cross(pos, FACE3_NEG)

    @assertEqual(19, pos % localID)

  end subroutine test_cross

  !!
  !! Test cell offset
  !!
@Test
  subroutine test_cellOffset()
    type(coord)                 :: pos
    real(defReal), dimension(3) :: ref
    real(defReal), parameter :: TOL = 1.0E-6_defReal

    ! ** 3D lattice, a selection of the 18 cells
    pos % localId = 1
    ref = [-3.0_defReal, -1.7320508_defReal, -2.5_defReal]
    @assertEqual(ref, uni1 % cellOffset(pos), TOL)

    pos % localId = 5
    ref = [ZERO, ZERO, -2.5_defReal]
    @assertEqual(ref, uni1 % cellOffset(pos), TOL)

    pos % localId = 18
    ref = [3.0_defReal, 1.7320508_defReal, 2.5_defReal]
    @assertEqual(ref, uni1 % cellOffset(pos), TOL)

    ! Background -> no offset
    pos % localId = 19
    ref = ZERO
    @assertEqual(ref, uni1 % cellOffset(pos), TOL)

    ! ** 2D lattice
    pos % localId = 2
    ref = [0.8660254_defReal, ZERO, -0.5_defReal]
    @assertEqual(ref, uni2 % cellOffset(pos), TOL)

    pos % localId = 5
    ref = ZERO
    @assertEqual(ref, uni2 % cellOffset(pos), TOL)

  end subroutine test_cellOffset

  !!
  !! Test getting a normal
  !!
@Test
  subroutine test_normal()
    type(coord)                 :: pos
    real(defReal), dimension(3) :: n
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    pos % r   = ZERO
    pos % dir = [ONE, ZERO, ZERO]

    n = uni1 % getNormal(FACE1_POS, pos)
    @assertEqual([0.5_defReal, 0.8660254_defReal, ZERO], n, TOL)

    n = uni1 % getNormal(FACE1_NEG, pos)
    @assertEqual([-0.5_defReal, -0.8660254_defReal, ZERO], n, TOL)

    n = uni1 % getNormal(FACE2_POS, pos)
    @assertEqual([-0.5_defReal, 0.8660254_defReal, ZERO], n, TOL)

    n = uni1 % getNormal(FACE2_NEG, pos)
    @assertEqual([0.5_defReal, -0.8660254_defReal, ZERO], n, TOL)

    n = uni1 % getNormal(FACE3_POS, pos)
    @assertEqual([-ONE, ZERO, ZERO], n, TOL)

    n = uni1 % getNormal(FACE3_NEG, pos)
    @assertEqual([ONE, ZERO, ZERO], n, TOL)

    n = uni1 % getNormal(AX_MIN, pos)
    @assertEqual([ZERO, ZERO, -ONE], n, TOL)

    n = uni1 % getNormal(AX_MAX, pos)
    @assertEqual([ZERO, ZERO, ONE], n, TOL)

    ! Outline surface: point on the +x face of the bounding box
    pos % r = [4.1547005373792452_defReal, ZERO, -2.5_defReal]
    pos % dir = [ONE, ZERO, ZERO]
    n = uni1 % getNormal(OUTLINE_SURF, pos)
    @assertEqual([ONE, ZERO, ZERO], n, TOL)

  end subroutine test_normal

end module hexUniverse_test