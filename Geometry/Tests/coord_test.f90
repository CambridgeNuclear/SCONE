module coord_test

  use numPrecision
  use coord_class, only : coord, coordList
  use funit

  implicit none

  ! Variables
  type(coordList) :: coords


contains

  !!
  !! Set up test envoroment
  !!
  !! Coords Placed at 3 levels
  !!   Level 1 -> global
  !!   Level 2 -> Rotated universe
  !!   Level 3 -> Translated universe
  !!
@Before
  subroutine set_up()
    real(defReal), dimension(3,3) :: mat

    ! Set Nesting
    coords % nesting = 3
    coords % matIdx = 2
    coords % uniqueId = 7

    ! Set Level 1
    coords % lvl(1) % r   = [1.0_defReal, 0.0_defReal, -1.0_defReal]
    coords % lvl(1) % dir = [ZERO, ONE, ZERO]
    coords % lvl(1) % uniIdx    = 1
    coords % lvl(1) % uniRootID = 1
    coords % lvl(1) % localID   = 1
    coords % lvl(1) % cellIdx   = 1

    ! Set Level 2
    ! Rotation Y -> Z; Z -> -Y
    mat = ZERO
    mat(1,1) = ONE
    mat(3,2) = -ONE
    mat(2,3) = ONE
    coords % lvl(2) % r   = [1.0_defReal, 0.0_defReal, -1.0_defReal]
    coords % lvl(2) % dir = [ZERO, ZERO, -ONE]
    coords % lvl(2) % uniIdx    = 2
    coords % lvl(2) % uniRootID = 6
    coords % lvl(2) % localID   = 3
    coords % lvl(2) % isRotated = .true.
    coords % lvl(2) % rotMat    = mat
    coords % lvl(2) % cellIdx   = 3

    ! Set Level 3
    ! Translation to origin
    coords % lvl(3) % r = ZERO
    coords % lvl(3) % dir = [ZERO, ZERO, -ONE]
    coords % lvl(3) % uniIdx    = 4
    coords % lvl(3) % uniRootID = 12
    coords % lvl(3) % localID   = 2
    coords % lvl(3) % cellIdx   = 0

  end subroutine set_up

  !!
  !! Clean test enviroment
  !!
@After
  subroutine clean_up()

    call coords % kill()

  end subroutine clean_up

  !!
  !! Test state changing procedures
  !!
@Test
  subroutine test_changing_state()

    ! Test State
    @assertTrue(coords % isPlaced())
    @assertFalse(coords % isAbove())
    @assertFalse(coords % isUninitialised())

    ! Change to above
    call coords % takeAboveGeom()
    @assertFalse(coords % isPlaced())
    @assertTrue(coords % isAbove())
    @assertFalse(coords % isUninitialised())

    ! Chenge to uninitialised
    call coords % kill()
    @assertFalse(coords % isPlaced())
    @assertFalse(coords % isAbove())
    @assertTrue(coords % isUninitialised())

  end subroutine test_changing_state

  !!
  !! Test nesting level changes & cell inquiry
  !!
@Test
  subroutine test_nesting_level()

    ! Move deeper
    call coords % addLevel()
    @assertEqual(4, coords % nesting)
    @assertEqual(0, coords % cell())

    ! Move to higher level
    call coords % decreaseLevel(2)
    @assertEqual(2, coords % nesting)
    @assertEqual(3, coords % cell())

  end subroutine test_nesting_level

  !!
  !! Test rotation
  !!
  !! Verifies only the deflection by mu !
  !!
@Test
  subroutine test_rotation()
    real(defReal)               :: mu, phi
    real(defReal), dimension(3) :: u1, u2, u3
    real(defReal), parameter    :: TOL = 1.0E-7_defReal

    mu = 0.3_defReal
    phi = 2.1_defReal

    ! Save pre-rotation direction
    u1 = coords % lvl(1) % dir
    u2 = coords % lvl(2) % dir
    u3 = coords % lvl(3) % dir

    call coords % rotate(mu, phi)

    ! Verify deflection
    @assertEqual(mu, dot_product(u1, coords % lvl(1) % dir), TOL)
    @assertEqual(mu, dot_product(u2, coords % lvl(2) % dir), TOL)
    @assertEqual(mu, dot_product(u3, coords % lvl(3) % dir), TOL)

  end subroutine test_rotation

  !!
  !! Test direction assigment
  !!
@Test
  subroutine test_direction_assigment()
    real(defReal), dimension(3) :: u1, u2, u3
    real(defReal), parameter    :: TOL = 1.0E-7_defReal

    ! Save pre-rotation direction
    u1 = coords % lvl(1) % dir
    u2 = coords % lvl(2) % dir
    u3 = coords % lvl(3) % dir

    ! Invert direction
    call coords % assignDirection(-coords % lvl(1) % dir)

    ! Verify
    @assertEqual(-u1, coords % lvl(1) % dir)
    @assertEqual(-u2, coords % lvl(2) % dir)
    @assertEqual(-u3, coords % lvl(3) % dir)

  end subroutine test_direction_assigment

  !!
  !! Test Movment
  !!
@Test
  subroutine test_movement()
    real(defReal), dimension(3) :: u1, u2, u3, r1, r2, r3
    real(defReal)               :: d
    real(defReal), parameter    :: TOL = 1.0E-7_defReal

    ! Move local
    d = 0.3_defReal
    r1 = coords % lvl(1) % r
    r2 = coords % lvl(2) % r
    r3 = coords % lvl(3) % r
    u1 = coords % lvl(1) % dir
    u2 = coords % lvl(2) % dir
    u3 = coords % lvl(3) % dir

    call coords % moveLocal(d, 3)

    ! Verify
    @assertEqual(r1 + d*u1, coords % lvl(1) % r, TOL)
    @assertEqual(r2 + d*u2, coords % lvl(2) % r, TOL)
    @assertEqual(r3 + d*u3, coords % lvl(3) % r, TOL)
    @assertTrue(coords % isPlaced())

    ! Move Global
    d = -13.0_defReal
    r1 = coords % lvl(1) % r
    u1 = coords % lvl(1) % dir

    call coords % moveGlobal(d)

    ! Verify
    @assertEqual(r1 + d*u1, coords % lvl(1) % r, TOL)
    @assertTrue(coords % isAbove())

  end subroutine test_movement

  !!
  !! Test coord validation
  !!
@Test
  subroutine test_coord_valid()

    @assertTrue(coords % lvl(1:coords % nesting) % isValid())

    call coords % kill()

    @assertFalse(coords % lvl % isValid())


  end subroutine test_coord_valid

end module coord_test
