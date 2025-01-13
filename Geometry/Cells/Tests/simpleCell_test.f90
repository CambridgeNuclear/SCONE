module simpleCell_test

  use numPrecision
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use surfaceShelf_class, only : surfaceShelf
  use simpleCell_class,   only : simpleCell
  use funit

  implicit none

  ! Parameters
  character(*), parameter :: SURF_DEF = "&
  & surf1 { id 13; type sphere; origin (0.0 0.0 0.0); radius 2.0;} &
  & surf2 { id 4;  type xPlane; x0 0.0;} &
  & surf3 { id 99; type yPlane; y0 0.0;}"

  ! Note that fill is not really needed to build a cell. It is used by cellShelf only
  ! Also note that the init procedure does not actually use the type, so it is left out.
  character(*), parameter :: CELL_DEF = "&
  & id 2; surfaces (-13 4 99 ); filltype outside; "


  ! Variables
  type(surfaceShelf) :: surfs
  type(simpleCell)   :: cell

contains

  !!
  !! Build the cell
  !!
@Before
  subroutine setUp()
    type(dictionary) :: dict

    call charToDict(dict, SURF_DEF)
    call surfs % init(dict)
    call dict % kill()
    call charToDict(dict, CELL_DEF)
    call cell % init(dict, surfs)

  end subroutine setUp

  !!
  !! Clean after tests
  !!
@After
  subroutine cleanUp()

    call surfs % kill()
    call cell % kill()

  end subroutine cleanUp

  !!
  !! Test Miscellaneous functionality
  !!
@Test
  subroutine test_misc()

    ! Test Id
    @assertEqual(2, cell % id())

    call cell % setID(7)
    @assertEqual(7, cell % id())

  end subroutine test_misc


  !!
  !! Test inside/outside determination
  !!
@Test
  subroutine test_inside()
    real(defReal), dimension(3) :: r, u

    ! Few points inside
    r = [1.0_defReal, 0.3_defReal, 0.1_defReal]
    u = [ONE, ZERO, ZERO]
    @assertTrue( cell % inside(r, u))

    r = [0.3_defReal, 1.3_defReal, 0.4_defReal]
    @assertTrue( cell % inside(r, u))

    ! Few points outside
    r = [-0.1_defReal, 1.2_defReal, 0.1_defReal]
    @assertFalse( cell % inside(r, u))

    r = [1.8_defReal, 0.8_defReal, 0.8_defReal]
    @assertFalse( cell % inside(r, u))

  end subroutine test_inside

  !!
  !! Test distance calculations
  !!
@Test
  subroutine test_distance()
    real(defReal), dimension(3) :: r, u
    real(defReal)               :: ref, d
    integer(shortInt)           :: idx, idx_ref
    real(defReal), parameter :: TOL = 1.0E-6

    ! Point inside
    ! X-Plane hit
    r = [0.3_defReal, 0.4_defReal, 0.0_defReal]
    u = [-ONE, ZERO, ZERO]
    ref = 0.3_defReal
    idx_ref = surfs % getIdx(4)
    call cell % distance(d, idx, r, u)
    @assertEqual(ref, d, TOL * ref)
    @assertEqual(idx_ref, idx)

    ! Y-Plane hit
    u = [ZERO, -ONE, ZERO]
    ref = 0.4_defReal
    idx_ref = surfs % getIdx(99)
    call cell % distance(d, idx, r, u)
    @assertEqual(ref, d, TOL * ref)
    @assertEqual(idx_ref, idx)

    ! Sphere hit
    u = [ONE, ZERO, ZERO]
    ref = 1.659591794_defReal
    idx_ref = surfs % getIdx(13)
    call cell % distance(d, idx, r, u)
    @assertEqual(ref, d, TOL * ref)
    @assertEqual(idx_ref, idx)

  end subroutine test_distance

end module simpleCell_test
