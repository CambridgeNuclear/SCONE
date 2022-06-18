module azimPinUniverse_test

  use numPrecision
  use universalVariables,    only : INF, SURF_TOL
  use dictionary_class,      only : dictionary
  use dictParser_func,       only : charToDict
  use charMap_class,         only : charMap
  use coord_class,           only : coord
  use surfaceShelf_class,    only : surfaceShelf
  use cellShelf_class,       only : cellShelf
  use azimPinUniverse_class, only : azimPinUniverse, MOVING_IN, MOVING_OUT, MOVING_CLOCK, &
                                    MOVING_ANTI, MOVING_CLOCK_FORWARD, MOVING_CLOCK_BACK
  use pfUnit_mod
  implicit none

  ! Parameters
  character(*), parameter :: UNI_DEF1 = &
  "id 7; type azimPinUniverse; naz 4; origin (0.0 0.0 0.0); rotation (0.0 0.0 0.0); &
  &radii (2.5 1.5 0.0); fills (u<7> u<14> void);"
  character(*), parameter :: UNI_DEF2 = &
  "id 8; type azimPinUniverse; naz 8; origin (0.0 0.0 0.0); rotation (0.0 0.0 0.0); &
  &radii (4.0 2.5 1.5 0.0); fills (u<20> u<7> u<14> void);"

  ! Variables
  type(surfaceShelf)    :: surfs
  type(cellShelf)       :: cells
  type(charMap)         :: mats
  type(azimPinUniverse) :: uni1, uni2


contains

  !!
  !! Set-up test enviroment
  !!
@Before
  subroutine setup()
    character(nameLen)                           :: name
    integer(shortInt), dimension(:), allocatable :: fill
    type(dictionary)                             :: dict
    integer(shortInt), dimension(:), allocatable :: fillArray

    ! Load void material
    name = 'void'
    call mats % add(name, 13)

    ! Build universe
    call charToDict(dict, UNI_DEF1)
    call uni1 % init(fill, dict, cells, surfs, mats)
    
    ! Set index
    call uni1 % setIdx(3)
    
    ! Verify fill array
    @assertEqual([-14, -14, -14, -14, -7, -7, -7, -7, 13, 13, 13, 13], fill)
    
    ! Build second universe
    fillArray = [-14, -14, -14, -14, -14, -14, -14, -14,-7, -7, -7, -7, -7, -7, -7, -7, &
            -20, -20, -20, -20, -20, -20, -20, -20, 13, 13, 13, 13, 13, 13, 13, 13]
    call charToDict(dict, UNI_DEF2)
    call uni2 % init(fill, dict, cells, surfs, mats)
    call uni2 % setIdx(26)
    @assertEqual(fillArray, fill)

  end subroutine setup

  !!
  !! Clean after test
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
  !! Test miscellaneous functionality
  !!
@Test
  subroutine test_misc()
    real(defReal), dimension(3,3) :: mat

    ! Get id
    @assertEqual(7, uni1 % id())

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
    real(defReal), dimension(3) :: r_ref, u_ref, r, dir
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    ! ** Enter into local cell 1
    r = [1.0_defReal, 0.0_defReal, 0.0_defReal ]
    dir = [ZERO, ZERO, ONE]

    call uni1 % enter(new, r, dir)

    ! Verify location
    r_ref = r
    u_ref = dir
    @assertEqual(r_ref, new % r, TOL)
    @assertEqual(u_ref, new % dir, TOL)
    @assertEqual(3, new % uniIdx)
    @assertEqual(1, new % localID)
    @assertEqual(0, new % cellIdx)
    
    ! ** Enter into local cell 2
    r = [0.0_defReal, 1.0_defReal, 0.0_defReal ]
    dir = [ZERO, ZERO, ONE]

    call uni1 % enter(new, r, dir)

    ! Verify location
    r_ref = r
    u_ref = dir
    @assertEqual(r_ref, new % r, TOL)
    @assertEqual(u_ref, new % dir, TOL)
    @assertEqual(3, new % uniIdx)
    @assertEqual(2, new % localID)
    @assertEqual(0, new % cellIdx)

    ! ** Enter into local cell 8
    r = [0.0_defReal, -2.3_defReal, -980.0_defReal ]
    dir = [ZERO, ZERO, ONE]

    call uni1 % enter(new, r, dir)

    ! Verify location
    r_ref = r
    u_ref = dir
    @assertEqual(r_ref, new % r, TOL)
    @assertEqual(u_ref, new % dir, TOL)
    @assertEqual(3, new % uniIdx)
    @assertEqual(8, new % localID)
    @assertEqual(0, new % cellIdx)

    ! ** Enter into local cell 11
    r = [-2.6_defReal, 0.0_defReal, -980.0_defReal ]
    dir = [ZERO, ZERO, ONE]

    call uni1 % enter(new, r, dir)

    ! Verify location
    r_ref = r
    u_ref = dir
    @assertEqual(r_ref, new % r, TOL)
    @assertEqual(u_ref, new % dir, TOL)
    @assertEqual(3, new % uniIdx)
    @assertEqual(11, new % localID)
    @assertEqual(0, new % cellIdx)

    ! VERIFY THAT ROTATION IS NOT SET (all angles were 0.0)
    @assertFalse(new % isRotated)

    ! In universe 2, enter local cell 21
    r = [-3.0_defReal, 0.1_defReal, 32.0_defReal]
    call uni2 % enter(new, r, dir)
    
    ! Verify location
    r_ref = r
    u_ref = dir
    @assertEqual(r_ref, new % r, TOL)
    @assertEqual(u_ref, new % dir, TOL)
    @assertEqual(26, new % uniIdx)
    @assertEqual(21, new % localID)
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
    real(defReal), parameter :: TOL = 1.0E-6_defReal

    ! ** In local cell 1 distance to radial boundary
    pos % r = [1.0_defReal, 0.0_defReal, 0.0_defReal]
    pos % dir = [ONE, ZERO, ZERO]
    pos % uniIdx  = 3
    pos % cellIdx = 0
    pos % localId = 1

    call uni1 % distance(d, surfIdx, pos)

    ref = 0.5_defReal
    @assertEqual(ref, d, ref * tol)
    @assertEqual(MOVING_OUT, surfIdx)

    ! ** In local cell 1 distance to anti-clockwise boundary
    pos % r = [1.0_defReal, 0.0_defReal, 0.0_defReal]
    pos % dir = [-SQRT2_2, SQRT2_2, ZERO]
    pos % uniIdx  = 3
    pos % cellIdx = 0
    pos % localId = 1

    call uni1 % distance(d, surfIdx, pos)

    ref = SQRT2_2
    @assertEqual(ref, d, ref * tol)
    @assertEqual(MOVING_ANTI, surfIdx)

    ! ** In local cell 3 distance to clockwise boundary
    pos % r = [-1.0_defReal, 0.0_defReal, 0.0_defReal]
    pos % dir = [SQRT2_2, SQRT2_2, ZERO]
    pos % uniIdx  = 3
    pos % cellIdx = 0
    pos % localId = 3

    call uni1 % distance(d, surfIdx, pos)

    ref = SQRT2_2
    @assertEqual(ref, d, ref * tol)
    @assertEqual(MOVING_CLOCK, surfIdx)

    ! ** In local cell 4 distance to clockwise boundary
    ! ** Moves back around!
    pos % r = [0.0_defReal, -1.0_defReal, 0.0_defReal]
    pos % dir = [SQRT2_2, SQRT2_2, ZERO]
    pos % uniIdx  = 3
    pos % cellIdx = 0
    pos % localId = 4

    call uni1 % distance(d, surfIdx, pos)

    ref = SQRT2_2
    @assertEqual(ref, d, ref * tol)
    @assertEqual(MOVING_CLOCK_FORWARD, surfIdx)

    ! ** In outermost cell moving away
    pos % r = [2.0_defReal, 1.6_defReal, 0.0_defReal]
    pos % dir = [ONE, ZERO, ZERO]
    pos % localId = 9

    call uni1 % distance(d, surfIdx, pos)
    @assertEqual(INF, d)
    ! Surface momento is undefined -> No crossing

    ! In ordinary cell in-between
    pos % r = [0.0_defReal, 1.6_defReal, 0.0_defReal]
    pos % dir = [ZERO, -ONE, ZERO]
    pos % localId = 5

    call uni1 % distance(d, surfIdx, pos)
    ref = 0.1_defReal
    @assertEqual(ref, d, ref * tol)
    @assertEqual(MOVING_IN, surfIdx)

    ! Universe 2


  end subroutine test_distance

  !!
  !! Test cell-to cell crossing
  !!
@Test
  subroutine test_cross()
    type(coord)       :: pos
    integer(shortInt) :: idx
    real(defReal) :: eps

    ! Cross from cell 2 to cell 6
    eps = HALF * SURF_TOL
    pos % r   = [0.0_defReal, 1.5_defReal-eps, 0.0_defReal]
    pos % dir = [ZERO, ONE, ZERO]
    pos % uniIdx = 8
    pos % cellIdx = 0
    pos % localId = 2

    idx = MOVING_OUT
    call uni1 % cross(pos, idx)

    @assertEqual(6, pos % localId)

    ! Cross from cell 6 to cell 2
    eps = HALF * SURF_TOL
    pos % r   = [0.0_defReal, 1.5_defReal+eps, 0.0_defReal]
    pos % dir = [ZERO, -ONE, ZERO]

    idx = MOVING_IN
    call uni1 % cross(pos, idx)

    @assertEqual(2, pos % localId)

    ! Cross from cell 1 to cell 4
    eps = HALF * SURF_TOL
    pos % r   = [SQRT2_2, -SQRT2_2, ZERO]
    pos % dir = [-SQRT2_2, -SQRT2_2, ZERO]
    pos % r   = pos % r - eps * pos % dir
    pos % localID = 1

    idx = MOVING_CLOCK_BACK
    call uni1 % cross(pos, idx)

    @assertEqual(4, pos % localId)

    ! Universe 2
    ! 

  end subroutine test_cross

  !!
  !! Test cell offset
  !!
@Test
  subroutine test_cellOffset()
    type(coord)       :: pos

    ! Cell 2
    pos % r   = [0.0_defReal, 1.0_defReal, 0.0_defReal]
    pos % dir = [ZERO, ONE, ZERO]
    pos % uniIdx = 3
    pos % cellIdx = 0
    pos % localId = 2

    @assertEqual([ZERO, ZERO, ZERO], uni1 % cellOffset(pos) )

    ! Cell 11
    pos % r   = [-7.0_defReal, 0.0_defReal, 0.0_defReal]
    pos % dir = [ZERO, ONE, ZERO]
    pos % uniIdx = 3
    pos % cellIdx = 0
    pos % localId = 11

    @assertEqual([ZERO, ZERO, ZERO], uni1 % cellOffset(pos) )

  end subroutine test_cellOffset

  !!
  !! Test surface transitions
  !!
  !! Check that there is no problem with distance calculations
  !! if particle is placed very close to an annulus or plane surface 
  !! (within SURF_TOL)
  !!
@Test
  subroutine test_edgeCases()
    type(coord)       :: pos
    integer(shortInt) :: idx, localID, cellIdx
    real(defReal)     :: eps, d
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    ! At boundary between cell 2 and 6
    eps = HALF * SURF_TOL
    pos % r   = [0.0_defReal, 1.5_defReal-eps, 0.0_defReal]
    pos % dir = [ONE, -0.00001_defReal, ZERO]
    pos % dir = pos % dir / norm2(pos % dir)
    pos % uniIdx = 8
    pos % cellIdx = 0

    ! Should find particle in cell 2
    ! And return very small distance -> MOVING OUT
    call uni1 % findCell(localID, cellIdx, pos % r, pos % dir)
    @assertEqual(2, localID)

    pos % localID = 2
    call uni1 % distance(d, idx, pos)

    @assertEqual(ZERO, d, 1.0E-3_defReal)
    @assertEqual(MOVING_OUT, idx)

    ! At boundary between cell 4 and 1
    eps = 1.1*SURF_TOL
    pos % r   = [SQRT2_2, -SQRT2_2, ZERO]
    pos % dir = [-SQRT2_2, -SQRT2_2, ZERO]
    pos % r   = pos % r + eps * pos % dir
    pos % dir = -pos % dir
    
    ! Should find particle in cell 4
    ! And return very small distance -> MOVING_CLOCK_FORWARD
    call uni1 % findCell(localID, cellIDx, pos % r, pos % dir)
    @assertEqual(4, localID)

    pos % localID = 4
    call uni1 % distance(d, idx, pos)

    @assertEqual(ZERO, d, 1.0E-3_defReal)
    @assertEqual(MOVING_CLOCK_FORWARD, idx)
    
  end subroutine test_edgeCases


end module azimPinUniverse_test
