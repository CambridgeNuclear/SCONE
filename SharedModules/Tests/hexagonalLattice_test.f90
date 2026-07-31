module hexagonalLattice_test

  use numPrecision
  use genericProcedures
  use dictionary_class,       only : dictionary
  use dictParser_func,        only : charToDict
  use hexagonalLattice_class
  use funit

  implicit none

  ! Parameters
  ! lat1: 3D lattice, "point" orientation, hexagonal plane in x-y, stacked along z
  !   halfwidth = 1.0 (pitch 2.0), axial pitch = 5.0, shape 3 x 3 x 2
  character(*), parameter :: LAT1_DEF = &
  "type zHexLattice; orientation 1; pitch (2.0 5.0); shape (3 3 2); "

  ! lat2: 2D lattice (infinite along the axis), "flat" orientation, hexagonal plane in x-z,
  !   stacked (infinitely) along y, offset from the global origin
  character(*), parameter :: LAT2_DEF = &
  "type yHexLattice; orientation 2; pitch (2.0 0.0); shape (2 2 0); origin (1.0 0.0 -1.0); "

  type(hexagonalLattice) :: lat1
  type(hexagonalLattice) :: lat2

contains

  !!
  !! Setup environment
  !!
@Before
  subroutine setUp()
    type(dictionary) :: dict

    call charToDict(dict, LAT1_DEF)
    call lat1 % init(dict)
    call dict % kill()

    call charToDict(dict, LAT2_DEF)
    call lat2 % init(dict)
    call dict % kill()

  end subroutine setUp

  !!
  !! Clean environment
  !!
@After
  subroutine clean()

    call lat1 % kill()
    call lat2 % kill()

  end subroutine clean

  !!
  !! Test basic getters
  !!
@Test
  subroutine test_getters()
    integer(shortInt), dimension(3) :: sizeN
    real(defReal), dimension(3) :: origin
    real(defReal), parameter :: TOL = 1.0E-9_defReal

    sizeN = lat1 % getSize()
    @assertEqual([3, 3, 2], sizeN)

    origin = lat1 % getOrigin()
    @assertEqual([ZERO, ZERO, ZERO], origin, TOL)

    @assertEqual(19, lat1 % getOutID())

    sizeN = lat2 % getSize()
    @assertEqual([2, 2, 1], sizeN)

    origin = lat2 % getOrigin()
    @assertEqual([1.0_defReal, ZERO, -1.0_defReal], origin, TOL)

    @assertEqual(5, lat2 % getOutID())

  end subroutine test_getters

  !!
  !! Test findCell
  !!
@Test
  subroutine test_findCell()
    real(defReal), dimension(3) :: r, u

    ! Centre of cell (2,2,1) -> localID 5
    r = [ZERO, ZERO, -2.5_defReal]
    u = [ZERO, ZERO, ONE]
    @assertEqual(5, lat1 % findCell(r, u))

    ! Well outside the lattice -> background
    r = [100.0_defReal, ZERO, -2.5_defReal]
    u = [ONE, ZERO, ZERO]
    @assertEqual(19, lat1 % findCell(r, u))

    ! lat2: centre of cell (2,1) -> localID 2 (offset by lat2's origin)
    r = [1.0_defReal + 0.8660254037844386_defReal, 123.4_defReal, -1.0_defReal - 0.5_defReal]
    u = [ZERO, ONE, ZERO]
    @assertEqual(2, lat2 % findCell(r, u))

  end subroutine test_findCell

  !!
  !! Test distance calculation
  !!
@Test
  subroutine test_distance()
    real(defReal), dimension(3) :: r, u
    real(defReal)     :: d, ref
    integer(shortInt) :: surfIdx
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    ! Well inside cell 5 (centre), moving along the axis -> hits the axial max face
    r = [ZERO, ZERO, -2.5_defReal]
    u = [ZERO, ZERO, ONE]
    call lat1 % distance(d, surfIdx, r, u)
    ref = 2.5_defReal
    @assertEqual(ref, d, TOL * ref)
    @assertEqual(AX_MAX, surfIdx)

    ! Well inside cell 5, moving toward the +x hexagon edge (face 3, negative side)
    r = [ZERO, ZERO, -2.5_defReal]
    u = [ONE, ZERO, ZERO]
    call lat1 % distance(d, surfIdx, r, u)
    ref = ONE
    @assertEqual(ref, d, TOL * ref)
    @assertEqual(FACE3_NEG, surfIdx)

    ! From well outside the lattice, aimed at it: the pure lattice component (unlike the
    ! hexLatUniverse wrapper, which uses a padded box for a single large jump) always
    ! resolves the nearest tiling cell directly and returns the small per-cell step - this
    ! is expected behaviour, not a bug: efficient long-range background stepping
    ! is the wrapper's responsibility.
    r = [100.0_defReal, ZERO, -2.5_defReal]
    u = [-ONE, ZERO, ZERO]
    call lat1 % distance(d, surfIdx, r, u)
    ref = ONE
    @assertEqual(ref, d, TOL * ref)
    @assertEqual(FACE3_POS, surfIdx)

  end subroutine test_distance

  !!
  !! Test cell offset
  !!
@Test
  subroutine test_getOffset()
    real(defReal), dimension(3) :: ref
    real(defReal), parameter :: TOL = 1.0E-6_defReal

    ref = [-3.0_defReal, -1.7320508_defReal, -2.5_defReal]
    @assertEqual(ref, lat1 % getOffset(1), TOL)

    ref = [ZERO, ZERO, -2.5_defReal]
    @assertEqual(ref, lat1 % getOffset(5), TOL)

    ! Background -> no offset
    ref = ZERO
    @assertEqual(ref, lat1 % getOffset(19), TOL)

    ! lat2 offset includes its own origin
    ref = [1.0_defReal + 0.8660254_defReal, ZERO, -1.5_defReal]
    @assertEqual(ref, lat2 % getOffset(2), TOL)

  end subroutine test_getOffset

  !!
  !! Test getting a normal
  !!
@Test
  subroutine test_getNormal()
    real(defReal), dimension(3) :: n
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    n = lat1 % getNormal(FACE1_POS)
    @assertEqual([0.5_defReal, 0.8660254_defReal, ZERO], n, TOL)

    n = lat1 % getNormal(FACE3_NEG)
    @assertEqual([ONE, ZERO, ZERO], n, TOL)

    n = lat1 % getNormal(AX_MIN)
    @assertEqual([ZERO, ZERO, -ONE], n, TOL)

    n = lat1 % getNormal(AX_MAX)
    @assertEqual([ZERO, ZERO, ONE], n, TOL)

  end subroutine test_getNormal

  !!
  !! Test self-consistency: findCell at the offset of every cell recovers that cell
  !!
@Test
  subroutine test_selfConsistency()
    integer(shortInt) :: id, foundID
    real(defReal), dimension(3) :: off, u

    u = [0.12345_defReal, 0.6789_defReal, 0.4567_defReal]
    u = u / norm2(u)

    do id = 1, 18
      off = lat1 % getOffset(id)
      foundID = lat1 % findCell(off, u)
      @assertEqual(id, foundID)
    end do

  end subroutine test_selfConsistency

end module hexagonalLattice_test
