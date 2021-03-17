module sphere_test

  use numPrecision
  use universalVariables, only : SURF_TOL, VACUUM_BC, REFLECTIVE_BC, INF
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use sphere_class,       only : sphere
  use pfUnit_mod

  implicit none

  character(*), parameter :: SPH_DEF = "id 7; origin (1.0 2.0 1.0); radius 2.0; "
  type(sphere)            :: surf

contains

  !!
  !! Build the surface
  !!
@Before
  subroutine setUp()
    type(dictionary) :: dict

    call charToDict(dict, SPH_DEF)
    call surf % init(dict)

  end subroutine setUp

  !!
  !! Clean after tests
  !!
@After
  subroutine cleanUp()

    call surf % kill()

  end subroutine cleanUp

  !!
  !! Test Misc functionality
  !!
@Test
  subroutine testMisc()
    real(defReal), dimension(6) :: aabb, ref
    real(defReal), parameter :: TOL = 1.0E-6_defReal

    ! Test ID
    @assertEqual(7, surf % id())

    ! Change ID
    call surf % setID(1)
    @assertEqual(1, surf % id())

    ! Bounding box
    ref = [-1.0_defReal, 0.0_defReal, -1.0_defReal, 3.0_defReal, 4.0_defReal, 3.0_defReal]
    aabb = surf % boundingBox()
    @assertEqual(ref, aabb, TOL)

    ! Name
    @assertEqual('sphere', surf % myType())

  end subroutine testMisc

  !!
  !! Test boundary conditions
  !!
@Test
  subroutine testBC()
    real(defReal), dimension(3) :: r, u, r_pre, u_pre

    ! Set Boundary conditions
    ! Should ignore extra entries
    call surf % setBC([VACUUM_BC, REFLECTIVE_BC, REFLECTIVE_BC])

    ! Apply BCs
    r = [1.0_defReal, 0.0_defReal, 1.0_defReal]
    u = [ZERO, ONE, ZERO]
    r_pre = r
    u_pre = u

    ! Explicit
    call surf % explicitBC(r, u)
    @assertEqual(r_pre, r)
    @assertEqual(u_pre, u)

    ! Transform
    call surf % transformBC(r, u)
    @assertEqual(r_pre, r)
    @assertEqual(u_pre, u)

  end subroutine testBC

  !!
  !! Test halfspace membership
  !!
@Test
  subroutine testHalfspace()
    real(defReal), dimension(3) :: r, u
    real(defReal)               :: eps

    ! Choose point At x-axis going at the surface
    r = [-1.0_defReal, 2.0_defReal, 1.0_defReal]
    u = [ONE, ZERO, ZERO]

    ! At the surface
    @assertFalse(surf % halfspace(r, u))

    ! Out within SURF_TOL
    eps = -0.5_defReal * SURF_TOL
    @assertFalse(surf % halfspace(r + eps*u, u))

    ! Out outside SURF_TOL
    eps = -1.00001_defReal * SURF_TOL
    @assertTrue(surf % halfspace(r + eps*u, u))

    ! Well Outside
    eps = -ONE
    @assertTrue(surf % halfspace(r + eps*u, u))

    ! Well within
    eps = ONE
    @assertFalse(surf % halfspace(r + eps*u, u))

    ! Tangent particle whould be outside
    u = [ZERO, ONE, ZERO]
    @assertTrue( surf % halfspace(r, u))

  end subroutine testHalfspace

  !!
  !! Test distance calculations
  !!
@Test
  subroutine testDistance()
    real(defReal), dimension(3) :: r, u
    real(defReal)               :: ref
    real(defReal), parameter :: TOL = 1.0E-7

    ! **Outside the sphere
    r = [-2.0_defReal, 2.0_defReal, 1.0_defReal]

    ! Perpendicular impact
    u = [ONE, ZERO, ZERO]
    ref = 1.0_defReal
    @assertEqual(ref, surf % distance(r, u), TOL * ref)

    ! Oblique impact
    ref = 1.5_defReal
    ! Calculate ange from Cosine theorem
    u(1) = (ref**2 + 3.0_defReal**2 - 2.0_defReal**2) / (TWO * ref * 3.0_defReal)
    u(2) = sqrt(ONE - u(1)**2)
    u = u / norm2(u)
    @assertEqual(ref, surf % distance(r, u), TOL * ref)

    ! **Exactly at the surface
    r  = [-1.0_defReal, 2.0_defReal, 1.0_defReal]

    ! Tangent particle
    u = [ZERO, ONE, ZERO]
    @assertEqual(INF, surf % distance(r, u))

    ! Particle going inside
    u = [ONE, ZERO, ZERO]
    ref = 4.0_defReal
    @assertEqual(ref, surf % distance(r, u))

    ! ** Outside in Surface tolerance
    r = r - [ONE, ZERO, ZERO] * HALF * SURF_TOL
    ref = 4.0_defReal + HALF * SURF_TOL
    @assertEqual(ref, surf % distance(r, u))

    ! **Inside the surface
    ! +ve direction
    r = [-0.5_defReal, 2.0_defReal, 1.0_defReal]
    u = [ONE, ZERO, ZERO]
    ref = 3.5_defReal
    @assertEqual(ref, surf % distance(r, u))

    ! -ve direction
    ref = 0.5_defReal
    @assertEqual(ref, surf % distance(r, -u))

  end subroutine testDistance

end module sphere_test
