module plane_test

  use numPrecision
  use universalVariables, only : SURF_TOL, VACUUM_BC, REFLECTIVE_BC, INF
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use plane_class,       only : plane
  use pfUnit_mod

  implicit none

  character(*), parameter :: PLANE_DEF = "id 7; coeffs (1.0 1.0 1.0 3.0);"
  type(plane)             :: surf

contains

  !!
  !! Build the surface
  !!
@Before
  subroutine setUp()
    type(dictionary) :: dict

    call charToDict(dict, PLANE_DEF)
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
    ref = [-INF, -INF, -INF, INF, INF, INF]
    aabb = surf % boundingBox()
    @assertEqual(ref, aabb, TOL)

    ! Name
    @assertEqual('plane', surf % myType())

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
    r = [1.0_defReal, 1.0_defReal, 1.0_defReal]
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
    real(defReal), dimension(3) :: r, u, u2
    real(defReal)               :: eps


    r = [ONE, ONE, ONE]
    u = [ONE, ONE, ONE]
    u = u /norm2(u)

    ! At the surface
    @assertTrue(surf % halfspace(r, u))

    ! In within SURF_TOL
    eps = -0.5_defReal * SURF_TOL * sqrt(3.0_defReal)
    @assertTrue(surf % halfspace(r + eps*u, u))

    ! In outside SURF_TOL
    eps = -1.00001_defReal * SURF_TOL * sqrt(3.0_defReal)
    @assertFalse(surf % halfspace(r + eps*u, u))

    ! Well Inside
    eps = -ONE
    @assertFalse(surf % halfspace(r + eps*u, u))

    ! Well Outside
    eps = ONE
    @assertTrue(surf % halfspace(r + eps*u, u))

    ! Tangent particle should use position
    u2 = [-ONE, ZERO, ONE]
    u2 = u2/norm2(u2)
    eps = HALF * SURF_TOL

    @assertTrue( surf % halfspace(r + eps*u, u2))
    @assertFalse( surf % halfspace(r - eps*u, u2))

  end subroutine testHalfspace

  !!
  !! Test distance calculations
  !!
@Test
  subroutine testDistance()
    real(defReal), dimension(3) :: r, u, u2
    real(defReal)               :: ref
    real(defReal), parameter :: TOL = 1.0E-7

    r = [-ONE, ZERO, ZERO]
    u = [ONE, ZERO, ZERO]

    ! ** Inside the surface
    ref = 4.0_defReal
    @assertEqual(ref, surf % distance(r, u), TOL * ref)

    ! ** Exactly at the surface
    r(1) = 3.0_defReal
    @assertEqual(INF, surf % distance(r, u))

    ! ** Inside within surface Tolerance
    r(1) = 3.0_defReal - SURF_TOL
    @assertEqual(INF, surf % distance(r, u))

    ! ** Inside close to the surface tolerance
    r(1) = 3.0_defReal - sqrt(3.0_defReal) * 1.1_defReal * SURF_TOL
    ref = sqrt(3.0_defReal) * 1.1_defReal * SURF_TOL
    ! Use very liberal tolerance becouse the distance is veary small and sensitive to FP precision
    @assertEqual(ref, surf % distance(r, u), 1.0E-2 * ref)

    ! ** Outside the surface
    r(1) = 4.0_defReal

    ! +ve Direction
    @assertEqual(INF, surf % distance(r, u))

    ! -ve Direction
    ref = 1.0_defReal
    @assertEqual(ref, surf % distance(r, -u), TOL * ref)

    ! ** Parallel to plane
    r = [-ONE, ZERO, ZERO]
    u2 = [-ONE, ZERO, ONE]
    u2 = u2/norm2(u2)
    @assertEqual(INF, surf % distance(r, u2))

  end subroutine testDistance

end module plane_test
