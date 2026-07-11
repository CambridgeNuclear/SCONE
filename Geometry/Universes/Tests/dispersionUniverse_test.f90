module dispersionUniverse_test

  use numPrecision
  use genericProcedures
  use universalVariables,       only : OVERLAP_MAT
  use dictionary_class,         only : dictionary
  use dictParser_func,          only : charToDict
  use charMap_class,            only : charMap
  use coord_class,              only : coord
  use surfaceShelf_class,       only : surfaceShelf
  use cellShelf_class,          only : cellShelf
  use cartesianLattice_class,   only : X_MIN, X_MAX, Y_MIN, Y_MAX, Z_MIN, Z_MAX
  use dispersionUniverse_class, only : dispersionUniverse, BACKGROUND_IDX, OUTSIDE_IDX, &
                                       OVERLAP_IDX, OUTLINE_SURF
  use funit

  implicit none

  ! Paths to the sphere-definition files consumed by readFile
  character(*), parameter :: FILE1 = '/tmp/dispersionUniverse_test_1.txt'
  character(*), parameter :: FILE2 = '/tmp/dispersionUniverse_test_2.txt'

  ! *** Universe 1: two well-separated, differently-sized spheres + background
  !     sphere 1 (fill=fuel):     radius 1, origin (-9, 2, 2)
  !     sphere 2 (fill=graphite): radius 5, origin (20,-2,-2)
  character(*), parameter :: UNI1_DEF = &
  "id 1; type dispersionUniverse; background void; checkOverlap 0; file "//FILE1//"; "

  ! *** Universe 2: two overlapping spheres, checkOverlap enabled
  !     sphere 1 (fill=fuel):     radius 2, origin (0, 0, 0)
  !     sphere 2 (fill=graphite): radius 2, origin (1, 0, 0)
  character(*), parameter :: UNI2_DEF = &
  "id 2; type dispersionUniverse; background void; checkOverlap 1; file "//FILE2//"; "

  ! Variables
  type(surfaceShelf)      :: surfs
  type(cellShelf)         :: cells
  type(charMap)           :: mats
  type(dispersionUniverse) :: uni1
  type(dispersionUniverse) :: uni2

contains

  !!
  !! Write the two sphere-definition files used by the universes above
  !!
  subroutine writeSphereFiles()
    integer(shortInt) :: unitNum

    open(newunit=unitNum, file=FILE1, status='replace', action='write')
    write(unitNum, '(4(F10.4,1X), A)')  1.0_defReal, -9.0_defReal,  2.0_defReal,  2.0_defReal, 'fuel'
    write(unitNum, '(4(F10.4,1X), A)')  5.0_defReal, 20.0_defReal, -2.0_defReal, -2.0_defReal, 'graphite'
    close(unitNum)

    open(newunit=unitNum, file=FILE2, status='replace', action='write')
    write(unitNum, '(4(F10.4,1X), A)')  2.0_defReal, 0.0_defReal, 0.0_defReal, 0.0_defReal, 'fuel'
    write(unitNum, '(4(F10.4,1X), A)')  2.0_defReal, 1.0_defReal, 0.0_defReal, 0.0_defReal, 'graphite'
    close(unitNum)

  end subroutine writeSphereFiles

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
    call mats % add(name, 2)
    name = 'fuel'
    call mats % add(name, 5)
    name = 'graphite'
    call mats % add(name, 8)

    call writeSphereFiles()

    ! Build universe 1
    call charToDict(dict, UNI1_DEF)
    call uni1 % init(fill, dict, cells, surfs, mats)
    call dict % kill()
    call uni1 % setIdx(4)

    ref = [5, 8, 2, 2, OVERLAP_MAT]
    @assertEqual(ref, fill)

    ! Build universe 2
    call charToDict(dict, UNI2_DEF)
    call uni2 % init(fill, dict, cells, surfs, mats)
    call dict % kill()
    call uni2 % setIdx(6)

    ref = [5, 8, 2, 2, OVERLAP_MAT]
    @assertEqual(ref, fill)

  end subroutine setUp

  !!
  !! Clean environment
  !!
@After
  subroutine clean()
    integer(shortInt) :: unitNum

    call surfs % kill()
    call cells % kill()
    call mats % kill()
    call uni1 % kill()
    call uni2 % kill()

    open(newunit=unitNum, file=FILE1, status='old', action='write')
    close(unitNum, status='delete')
    open(newunit=unitNum, file=FILE2, status='old', action='write')
    close(unitNum, status='delete')

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
  !! Test entering a universe
  !!
@Test
  subroutine test_enter()
    type(coord) :: new
    real(defReal), dimension(3) :: r, dir
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    ! Inside sphere 1 (centre)
    r   = [-9.0_defReal, 2.0_defReal, 2.0_defReal]
    dir = [ZERO, ZERO, ONE]

    call uni1 % enter(new, r, dir)

    @assertEqual(4, new % uniIdx)
    @assertEqual(1, new % localID)
    @assertEqual(0, new % cellIdx)

    ! Inside sphere 2 (centre)
    r   = [20.0_defReal, -2.0_defReal, -2.0_defReal]
    dir = [ONE, ZERO, ZERO]

    call uni1 % enter(new, r, dir)

    @assertEqual(2, new % localID)

    ! In the background, well inside the outline, away from both spheres
    r   = [ZERO, ZERO, ZERO]
    dir = [ONE, ZERO, ZERO]

    call uni1 % enter(new, r, dir)

    @assertEqual(2 + BACKGROUND_IDX, new % localID)

    ! Completely outside the outline box
    r   = [1000.0_defReal, 1000.0_defReal, 1000.0_defReal]
    dir = [ONE, ZERO, ZERO]

    call uni1 % enter(new, r, dir)

    @assertEqual(2 + OUTSIDE_IDX, new % localID)

    ! *** Universe 2: overlap detection
    ! In the region shared by both spheres
    r   = [0.5_defReal, ZERO, ZERO]
    dir = [ZERO, ZERO, ONE]

    call uni2 % enter(new, r, dir)

    @assertEqual(6, new % uniIdx)
    ! nSpheres + 2 -> OVERLAP_MAT cell
    @assertEqual(2 + OVERLAP_IDX, new % localID)

    ! Only inside sphere 1 (not sphere 2)
    r   = [-1.5_defReal, ZERO, ZERO]
    dir = [ZERO, ZERO, ONE]

    call uni2 % enter(new, r, dir)

    @assertEqual(1, new % localID)

  end subroutine test_enter

  !!
  !! Test distance calculation
  !!
@Test
  subroutine test_distance()
    real(defReal)     :: d, ref
    integer(shortInt) :: surfIdx
    type(coord)       :: pos
    real(defReal), dimension(3) :: dir
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    ! Inside sphere 1, starting exactly at its centre.
    ! Distance to the sphere surface from the centre is exactly the radius,
    ! regardless of direction or how the acceleration grid subdivides this region.
    dir = [ONE, ONE, ONE]
    dir = dir / norm2(dir)

    pos % r      = [-9.0_defReal, 2.0_defReal, 2.0_defReal]
    pos % dir    = dir
    pos % uniIdx = 4
    pos % cellIdx = 0
    pos % localId = 1

    call uni1 % distance(d, surfIdx, pos)

    ref = ONE
    @assertEqual(ref, d, TOL * ref)
    ! -s + OUTLINE_SURF, s = 1
    @assertEqual(-8, surfIdx)

    ! Point (2,0,0) sits in a grid cell containing no spheres at all
    ! (cell x-range ~[1.666, 7.500], well clear of both spheres).
    ! The +x wall of that cell sits at exactly x = 7.5 (by symmetry of the
    ! grid about the outline's centre), so the distance is exactly 5.5.
    pos % r      = [2.0_defReal, ZERO, ZERO]
    pos % dir    = [ONE, ZERO, ZERO]
    pos % uniIdx = 4
    pos % cellIdx = 0
    pos % localId = 3   ! background

    call uni1 % distance(d, surfIdx, pos)

    ref = 5.5_defReal
    @assertEqual(ref, d, TOL * ref)
    @assertEqual(X_MAX, surfIdx)

    ! Sphere 2 (radius 5, centre (20,-2,-2)) spans grid cells 5 and 6 in x.
    ! Start inside sphere 2, within cell 5 (x in [13.334, 19.168]):
    !   r = (16,-1,-1); distance from centre = sqrt(16+1+1) = sqrt(18) < 5, so inside.
    ! Moving +x, the sphere's true exit is at d = 4 + sqrt(23) =~ 8.796,
    ! but cell 5's +x wall sits at x =~ 19.16783, giving d =~ 3.16783 --
    ! strictly closer than the sphere exit, so the grid wall must win.
    pos % r      = [16.0_defReal, -1.0_defReal, -1.0_defReal]
    pos % dir    = [ONE, ZERO, ZERO]
    pos % uniIdx = 4
    pos % cellIdx = 0
    pos % localId = 2   ! inside sphere 2

    call uni1 % distance(d, surfIdx, pos)

    ref = 3.167833333333333_defReal
    @assertEqual(ref, d, TOL * ref)
    @assertEqual(X_MAX, surfIdx)

  end subroutine test_distance

  !!
  !! Test cell-to-cell crossing
  !!
@Test
  subroutine test_cross()
    type(coord)   :: pos
    real(defReal) :: eps

    eps = HALF * SURF_TOL

    ! Cross outward, out of sphere 1 into the background
    pos % r      = [-8.0_defReal + eps, 2.0_defReal, 2.0_defReal]
    pos % dir    = [ONE, ZERO, ZERO]
    pos % uniIdx = 4
    pos % cellIdx = 0
    pos % localId = 1

    call uni1 % cross(pos, -2)

    @assertEqual(3, pos % localID)

    ! Cross inward, from the background into sphere 1
    pos % r      = [-10.0_defReal + eps, 2.0_defReal, 2.0_defReal]
    pos % dir    = [ONE, ZERO, ZERO]
    pos % localId = 3

    call uni1 % cross(pos, -2)

    @assertEqual(1, pos % localID)

  end subroutine test_cross

  !!
  !! Test cell offset - should translate to the sphere
  !!
@Test
  subroutine test_cellOffset()
    type(coord)                 :: pos
    real(defReal), dimension(3) :: ref
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    pos % r      = [-9.0_defReal, 2.0_defReal, 2.0_defReal]
    pos % dir    = [ONE, ZERO, ZERO]
    pos % uniIdx = 4
    pos % cellIdx = 0
    pos % localId = 1

    ref = pos % r
    @assertEqual(ref, uni1 % cellOffset(pos), TOL)

  end subroutine test_cellOffset

  !!
  !! Test getting a normal
  !!
@Test
  subroutine test_normal()
    type(coord)                 :: pos
    real(defReal), dimension(3) :: n

    ! Normal on the outline's +x face
    pos % r   = [25.00175_defReal, -2.0_defReal, -2.0_defReal]
    pos % dir = [ONE, ZERO, ZERO]

    n = uni1 % getNormal(OUTLINE_SURF, pos)
    @assertEqual([1.0_defReal, 0.0_defReal, 0.0_defReal], n, 1.0E-7_defReal)

    ! Normal on sphere 1's surface (+x point) -> outward radial direction
    pos % r   = [-8.0_defReal, 2.0_defReal, 2.0_defReal]
    pos % dir = [ONE, ZERO, ZERO]

    n = uni1 % getNormal(-2, pos)
    @assertEqual([1.0_defReal, 0.0_defReal, 0.0_defReal], n, 1.0E-7_defReal)

  end subroutine test_normal

end module dispersionUniverse_test