module quadric_test

  use numPrecision
  use universalVariables, only : SURF_TOL, VACUUM_BC, REFLECTIVE_BC, INF
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use quadric_class,      only : quadric
  use funit

  implicit none

  ! z-cylinder with radius 2
  character(*), parameter :: CYL_DEF = "id 7; coeffs (1.0 1.0 0 0 0 0 0 0 0 -4); "
  ! Ellipse shifted by +1 in x
  character(*), parameter :: ELP_DEF = "id 3; coeffs (0.25 0.0625 0.0 0 0 0 -0.5 0 0 -0.75); "
  type(quadric)           :: cyl, elp

contains

  !!
  !! Build the surface
  !!
@Before
  subroutine setUp()
    type(dictionary) :: dict

    call charToDict(dict, CYL_DEF)
    call cyl % init(dict)
    call charToDict(dict, ELP_DEF)
    call elp % init(dict)

  end subroutine setUp

  !!
  !! Clean after tests
  !!
@After
  subroutine cleanUp()

    call cyl % kill()
    call elp % kill()

  end subroutine cleanUp

  !!
  !! Test Misc functionality
  !!
@Test
  subroutine testMisc()
    real(defReal), dimension(6) :: aabb, ref
    real(defReal), parameter :: TOL = 1.0E-6_defReal

    ! Test ID
    @assertEqual(7, cyl % id())

    ! Change ID
    call cyl % setID(1)
    @assertEqual(1, cyl % id())
    
    ! Bounding box
    ref = [-INF, -INF, -INF, INF, INF, INF]
    aabb = cyl % boundingBox()
    @assertEqual(ref, aabb, TOL)

    ! Name
    @assertEqual('quadric', cyl % myType())

  end subroutine testMisc

  !!
  !! Test halfspace membership
  !!
@Test
  subroutine testHalfspace()
    real(defReal), dimension(3) :: r, u
    real(defReal)               :: eps

    ! Cylinder tests
    ! Point inside
    r = [-0.5_defReal, 1.0_defReal, 10.0_defReal]
    u = [-ONE, ZERO, ZERO]
    @assertFalse(cyl % halfspace(r, u))

    ! At the surface
    r = [TWO, ZERO, 10.0_defReal]
    @assertFalse(cyl % halfspace(r, u))

    ! Out within SURF_TOL
    eps = -0.5_defReal * SURF_TOL
    @assertFalse(cyl % halfspace(r + eps*u, u))

    ! Out outside SURF_TOL
    eps = -8.0001 * SURF_TOL
    @assertTrue(cyl % halfspace(r + eps*u, u))

    ! Well Outside
    eps = -ONE
    @assertTrue(cyl % halfspace(r + eps*u, u))

    ! Well within
    eps = ONE
    @assertFalse(cyl % halfspace(r + eps*u, u))

    ! Tangent particle would be outside
    u = [ZERO, ONE, ZERO]
    @assertTrue( cyl % halfspace(r, u))

    ! Ellipse tests
    ! Centred on 1, minor radius of 2 (x),
    ! major radius of 4 (y)
    r = [ZERO, 2.0_defReal, -8.0_defReal]
    u = [-ONE, ZERO, ZERO]
    @assertFalse(elp % halfspace(r, u))
    
    r = [-ONE, 2.0_defReal, -8.0_defReal]
    @assertTrue(elp % halfspace(r, u))
    
    ! At the surface
    r = [3.0_defReal, ZERO, -8.0_defReal]
    @assertFalse(elp % halfspace(r, u))
    
    ! Out within SURF_TOL
    eps = -0.5_defReal * SURF_TOL
    @assertFalse(elp % halfspace(r + eps*u, u))
    
    ! Out outside SURF_TOL
    eps = -1.501_defReal * SURF_TOL
    @assertTrue(elp % halfspace(r + eps*u, u))
    
    ! Tangent particle would be outside
    u = [ZERO, ONE, ZERO]
    @assertTrue(elp % halfspace(r, u))

  end subroutine testHalfspace

  !!
  !! Test distance calculations
  !!
@Test
  subroutine testDistance()
    real(defReal), dimension(3) :: r, u
    real(defReal)               :: ref
    real(defReal), parameter :: TOL = 1.0E-7

    ! **Outside the cylinder
    r = [-4.0_defReal, 0.0_defReal, 1.0_defReal]

    ! Perpendicular impact
    u = [ONE, ZERO, ZERO]
    ref = 2.0_defReal
    @assertEqual(ref, cyl % distance(r, u), TOL * ref)

    ! Oblique impact
    r = [-3.0_defReal, 0.5_defReal, 30.0_defReal]
    u = [0.744208408_defReal, 0.248069469_defReal, 0.620173673_defReal]
    ref = 1.634986685_defReal
    @assertEqual(ref, cyl % distance(r, u), TOL * ref)

    ! **Exactly at the surface
    r  = [0.0_defReal, 2.0_defReal, 1.0_defReal]

    ! Tangent particle
    u = [ONE, ZERO, ZERO]
    @assertEqual(INF, cyl % distance(r, u))

    ! Particle inside going out
    u = [0.707106781187_defReal, 0.707106781187_defReal, ZERO]
    r  = [-0.5_defReal, -1.0_defReal, 1.0_defReal]
    ref = 3.0291621403_defReal
    @assertEqual(ref, cyl % distance(r, u), TOL * ref)

    ! ** Outside in Surface tolerance
    r  = [-2.0_defReal, 0.0_defReal, 1.0_defReal]
    u = [ONE, ZERO, ZERO]
    r = r - u * HALF * SURF_TOL
    ref = 4.0_defReal + HALF * SURF_TOL
    @assertEqual(ref, cyl % distance(r, u), TOL * ref)

    ! And now for the ellipse
    ! Outside going in
    r = [0.0_defReal, -6.0_defReal, 12.0_defReal]
    u = [ZERO, ONE, ZERO]
    ref = 2.53589838486_defReal
    @assertEqual(ref, elp % distance(r, u), TOL * ref)
    
    ! Same but a less easy angle
    r = [-1.0_defReal, -4.0_defReal, 12.0_defReal]
    u = [0.707106781187_defReal, 0.707106781187_defReal, ZERO]
    ref = 1.13137084990_defReal
    @assertEqual(ref, elp % distance(r, u), TOL * ref)

    ! Inside going out
    r = [0.0_defReal, -2.0_defReal, 12.0_defReal]
    ref = 4.16282187604_defReal
    @assertEqual(ref, elp % distance(r, u), TOL * ref)

  end subroutine testDistance

end module quadric_test
