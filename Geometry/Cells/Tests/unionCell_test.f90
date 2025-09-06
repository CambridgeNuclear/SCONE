module unionCell_test

  use numPrecision
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use surfaceShelf_class, only : surfaceShelf
  use unionCell_class,    only : unionCell
  use funit

  implicit none

  ! Parameters
  character(*), parameter :: SURF_DEF = "&
  & surf1 { id 13; type sphere; origin (0.0 0.0 0.0); radius 5.0;} &
  & surf2 { id 15; type box; origin (0.0 0.0 0.0); halfwidth ( 4.6 4.6 4.6);} &
  & surf3 { id 4;  type xTruncCylinder; origin (0 0 0); halfwidth 5.1; radius 1.0;} &
  & surf4 { id 5;  type yTruncCylinder; origin (0 0 0); halfwidth 5.1; radius 1.0;} &
  & surf5 { id 6;  type zTruncCylinder; origin (0 0 0); halfwidth 5.1; radius 1.0;} "

  ! Note that fill is not really needed to build a cell. It is used by cellShelf only
  ! Also note that the init procedure does not actually use the type, so it is left out.
  character(*), parameter :: EASYCELL_DEF = "&
  & id 2; surfaces [-13 -15 ]; filltype outside; "
  character(*), parameter :: COMPLEXCELL_DEF = "&
  & id 3; surfaces [< -13 -15 > # < -4 : -5 : -6 > ]; filltype outside; "


  ! Variables
  type(surfaceShelf) :: surfs
  type(unionCell)    :: easyCell, complexCell

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

    call charToDict(dict, EASYCELL_DEF)
    call easyCell % init(dict, surfs)
    call dict % kill()

    call charToDict(dict, COMPLEXCELL_DEF)
    call complexCell % init(dict, surfs)

  end subroutine setUp

  !!
  !! Clean after tests
  !!
@After
  subroutine cleanUp()

    call surfs % kill()
    call easyCell % kill()
    call complexCell % kill()

  end subroutine cleanUp

  !!
  !! Test Miscellaneous functionality
  !!
@Test
  subroutine test_misc()

    ! Test Id
    @assertEqual(3, complexCell % id())

    call complexCell % setID(7)
    @assertEqual(7, complexCell % id())

  end subroutine test_misc


  !!
  !! Test inside/outside determination
  !!
@Test
  subroutine test_inside()
    real(defReal), dimension(3) :: r, u

    u = [ONE, ZERO, ZERO]
    
    ! Few points inside
    r = [-4.59_defReal, 0.0_defReal, 0.0_defReal]
    @assertTrue( easyCell % inside(r, u))

    r = [0.3_defReal, 1.3_defReal, 0.4_defReal]
    @assertTrue( easyCell % inside(r, u))

    r = [3.0_defReal, 2.3_defReal, 0.4_defReal]
    @assertTrue( complexCell % inside(r, u))

    r = [-1.1_defReal, 1.05_defReal, 0.0_defReal]
    @assertTrue( complexCell % inside(r, u))

    ! Few points outside
    r = [0.0_defReal, 4.9_defReal, 0.0_defReal]
    @assertFalse( easyCell % inside(r, u))

    r = [4.58_defReal, 4.58_defReal, -4.58_defReal]
    @assertFalse( easyCell % inside(r, u))
    
    r = [0.1_defReal, 0.3_defReal, -0.2_defReal]
    @assertFalse( complexCell % inside(r, u))
    
    r = [-4.59_defReal, 0.0_defReal, 0.25_defReal]
    @assertFalse( complexCell % inside(r, u))
    
    r = [4.58_defReal, 4.58_defReal, -4.58_defReal]
    @assertFalse( complexCell % inside(r, u))

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

    ! Points inside
    r = [0.2_defReal, -0.1_defReal, 0.0_defReal]
    u = [-0.980580676_defReal, ZERO, 0.196116135_defReal]
    ref = 4.895058733_defReal
    idx_ref = surfs % getIdx(15)
    call easyCell % distance(d, idx, r, u)
    @assertEqual(ref, d, TOL * ref)
    @assertEqual(idx_ref, idx)

    ! Point outside - should hit an imaginary sphere surface
    r = [6.0_defReal, 0.0_defReal, 0.0_defReal]
    u = [-ONE, ZERO, ZERO]
    ref = ONE
    idx_ref = surfs % getIdx(13)
    call easyCell % distance(d, idx, r, u)
    @assertEqual(ref, d, TOL * ref)
    @assertEqual(idx_ref, idx)
    
    ! Point outside - should hit a real sphere surface
    r = [4.5_defReal, 4.5_defReal, 4.5_defReal]
    u = -0.577350269_defReal
    ref = 2.794228634_defReal
    idx_ref = surfs % getIdx(13)
    call easyCell % distance(d, idx, r, u)
    @assertEqual(ref, d, TOL * ref)
    @assertEqual(idx_ref, idx)

    ! Cylinder hit
    r = [1.2_defReal, 0.0_defReal, 0.0_defReal]
    u = [ZERO, -ONE, ZERO]
    ref = ONE
    idx_ref = surfs % getIdx(4)
    call complexCell % distance(d, idx, r, u)
    @assertEqual(ref, d, TOL * ref)
    @assertEqual(idx_ref, idx)

  end subroutine test_distance

end module unionCell_test
