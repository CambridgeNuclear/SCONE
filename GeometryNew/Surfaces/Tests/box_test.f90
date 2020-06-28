module box_test

  use numPrecision
  use universalVariables
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use box_class,          only : box
  use pfUnit_mod

  implicit none

  character(*), parameter :: BOX_DEF = "id 7; origin (1.0 2.0 1.0); halfwidth (1.0 2.0 3.0); "
  type(box)               :: surf

contains

  !!
  !! Build the surface
  !!
@Before
  subroutine setUp()
    type(dictionary) :: dict

    call charToDict(dict, BOX_DEF)
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
    ref = [0.0_defReal, 0.0_defReal, -2.0_defReal, 2.0_defReal, 4.0_defReal, 4.0_defReal]
    aabb = surf % boundingBox()
    @assertEqual(ref, aabb, TOL)

    ! Name
    @assertEqual('box', surf % myType())

  end subroutine testMisc

  !!
  !! Test boundary conditions
  !!
@Test
  subroutine testBC()
    real(defReal), dimension(3) :: r, u, r_ref, u_ref
    real(defReal), parameter :: TOL = 1.0E-6

    ! Assign boundary conditions
    call surf % setBC([VACUUM_BC, REFLECTIVE_BC, PERIODIC_BC, &
                      PERIODIC_BC, REFLECTIVE_BC, REFLECTIVE_BC])

    ! ** Explicit BC
    ! Vacuum surface
    r = [0.0_defReal, 1.0_defReal, -1.0_defReal]
    u = [ONE, ZERO, ZERO]
    r_ref = r
    u_ref = u
    call surf % explicitBC(r, u)
    @assertEqual(r_ref, r, TOL)
    @assertEqual(u_ref, u, TOL)

    ! Reflective surface
    r = [1.0_defReal, 1.0_defReal, -2.0_defReal]
    u = [ZERO, ZERO, -ONE]
    r_ref = r
    u_ref = -u
    call surf % explicitBC(r, u)
    @assertEqual(r_ref, r, TOL)
    @assertEqual(u_ref, u, TOL)

    ! Periodic surface
    r = [1.0_defReal, 0.0_defReal, -1.0_defReal]
    u = [ZERO, -ONE, ZERO]
    r_ref = [1.0_defReal, 4.0_defReal, -1.0_defReal]
    u_ref = u
    call surf % explicitBC(r, u)
    @assertEqual(r_ref, r, TOL)
    @assertEqual(u_ref, u, TOL)

    ! Try the corner
    r = [2.0_defReal, 0.0_defReal, 4.0_defReal]
    u = [ONE, -ONE, ONE]
    u = u / norm2(u)

    r_ref = [2.0_defReal, 4.0_defReal, 4.0_defReal]
    u_ref = [-ONE, -ONE, -ONE]
    u_ref = u_ref / norm2(u_ref)

    call surf % explicitBC(r, u)
    @assertEqual(r_ref, r, TOL)
    @assertEqual(u_ref, u, TOL)

    ! ** Transform BC
    !
    r = [4.5_defReal, -10.0_defReal, 4.0_defReal]
    u = [ONE, ONE, ONE]
    u = u / norm2(u)

    r_ref = [-0.5_defReal, 2.0_defReal, 4.0_defReal]
    u_ref = [-ONE, ONE, -ONE]
    u_ref = u_ref / norm2(u_ref)

    call surf % transformBC(r, u)
    @assertEqual(r_ref, r, TOL)
    @assertEqual(u_ref, u, TOL)

    ! Try the corner
    r = [2.0_defReal, 0.0_defReal, 4.0_defReal]
    u = [ONE, -ONE, ONE]
    u = u / norm2(u)

    r_ref = [2.0_defReal, 4.0_defReal, 4.0_defReal]
    u_ref = [-ONE, -ONE, -ONE]
    u_ref = u_ref / norm2(u_ref)

    call surf % transformBC(r, u)
    @assertEqual(r_ref, r, TOL)
    @assertEqual(u_ref, u, TOL)

  end subroutine testBC

  !!
  !! Test halfspace membership
  !!
@Test
  subroutine testHalfspace()
    real(defReal), dimension(3) :: r, u, u2
    real(defReal)               :: eps

    ! ** Well inside the box
    r = [0.5_defReal, 1.0_defReal, 3.6_defReal]
    u = [ZERO, ZERO, ONE]
    @assertFalse(surf % halfspace(r, u))

    ! Difrent octant
    r = [1.5_defReal, 1.0_defReal, 0.5_defReal]
    @assertFalse(surf % halfspace(r, u))

    ! ** Well outside the box
    ! Make sure point is at one of the planes outside
    r = [-0.5_defReal, 0.0_defReal, 3.6_defReal]
    @assertTrue(surf % halfspace(r, u))

    ! Diffrent octant
    r = [2.0_defReal, 5.0_defReal, 0.5_defReal]
    @assertTrue(surf % halfspace(r, u))

    ! ** Proximity of a surface
    r = [1.5_defReal, 4.0_defReal, 0.5_defReal]
    u = [ZERO, ONE, ZERO]
    @assertTrue(surf % halfspace(r, u))

    ! Inside, within surface tolerance
    eps = -HALF * SURF_TOL
    @assertTrue(surf % halfspace(r + eps*u, u))

    ! Inside, a bit outside surface tolerance
    eps = -1.00001_defReal * SURF_TOL
    @assertFalse(surf % halfspace(r + eps*u, u))

    ! Parallel to the surface
    u2 = [ONE, ZERO, ZERO]
    eps = HALF * SURF_TOL

    @assertTrue(surf % halfspace(r + eps*u, u2))
    @assertFalse(surf % halfspace(r - eps*u, u2))

  end subroutine testHalfspace

  !!
  !! Test distance calculations
  !!
@Test
  subroutine testDistance()
    real(defReal), dimension(3) :: r, u
    real(defReal)               :: ref
    real(defReal), parameter :: TOL = 1.0E-7

    ! ** Outside the box
    ! Note that point is in x-z plane (same as bottom of the box)
    r = [-2.0_defReal, 0.000_defReal, 1.0_defReal]

    ! Direct impact
    u = [ONE, ZERO, ZERO]
    ref = 2.0_defReal
    @assertEqual(ref, surf % distance(r, u), TOL*ref)

    ! Miss
    u = [ONE, -0.0001_defReal, ONE]
    u = u /norm2(u)
    @assertEqual(INF, surf % distance(r, u))

    ! Oblique Hit
    u = [ONE, ZERO, ONE]
    u = u /norm2(u)
    ref = 2.0_defReal * SQRT2
    @assertEqual(ref, surf % distance(r, u), TOL*ref)

    ! Corner skim
    ! Use dirty values
    r = [-2.0_defReal, 0.000_defReal, 1.0_defReal/3.0_defReal]
    u = [TWO, ZERO, 3.0_defReal + 2.0_defReal/3.0_defReal]
    u = u /norm2(u)
    @assertEqual(INF, surf % distance(r, u), TOL*ref)

    ! ** At the surface
    r = [ONE, TWO, -2.0_defReal]
    u = [ZERO, ZERO, ONE]
    ref = 6.0_defReal
    @assertEqual(ref, surf % distance(r, u), TOL*ref)
    @assertEqual(INF, surf % distance(r, -u))

    ! Outside within surface tolerance
    r = [ONE, TWO, -2.0_defReal - HALF * SURF_TOL]
    ref = 6.0_defReal + HALF * SURF_TOL
    @assertEqual(ref, surf % distance(r, u), TOL*ref)

    ! ** Inside
    r = [ONE, TWO, -1.0_defReal]
    u = [ZERO, ONE, ONE]
    u = u/norm2(u)

    ref = TWO * SQRT2
    @assertEqual(ref, surf % distance(r, u), TOL*ref)

    ref = SQRT2
    @assertEqual(ref, surf % distance(r, -u), TOL*ref)

  end subroutine testDistance

  !!
  !! Test Edge Cases
  !!
  !! Test unlikley cases to make sure that halfspace + distance
  !! procedures allow particle to escape
  !!
@Test
  subroutine testEdgeCases()
    real(defReal), dimension(3) :: r, u, u2
    real(defReal)               :: eps, d
    logical(defBool)            :: hs

    ! ** Corner
    ! * Particle is almost at the corner
    !   Either it is outside or can escape with a short movment in next step
    eps =  5.0_defReal * epsilon(eps)
    r = [2.0_defReal-eps, 0.0_defReal, 4.0_defReal-eps]
    u = [HALF, ZERO, -ONE]
    u = u/norm2(u)
    hs = surf % halfspace(r, u)
    if (.not.hs) then ! Perform small movment
      d = surf % distance(r, u)
      @assertTrue( abs(d) < 1.0E-6)
      @assertTrue(surf % halfspace(r + d*u, u))
    end if

    ! Try in other direction
    u = [-ONE, ZERO, HALF]
    u = u/norm2(u)
    hs = surf % halfspace(r, u)
    if (.not.hs) then ! Perform small movment
      d = surf % distance(r, u)
      @assertTrue( abs(d) < 1.0E-6)
      @assertTrue(surf % halfspace(r + d*u, u))
    end if

    ! Try asymertic corner
    r = [2.0_defReal-TWO*eps, 0.0_defReal, 4.0_defReal-eps]
    u = [HALF, ZERO, -ONE]
    u = u/norm2(u)
    hs = surf % halfspace(r, u)
    if (.not.hs) then ! Perform small movment
      d = surf % distance(r, u)
      @assertTrue( abs(d) < 1.0E-6)
      @assertTrue(surf % halfspace(r + d*u, u))
    end if

    ! Asymetric corner position
    ! Try other direction
    r = [2.0_defReal-eps, 0.0_defReal, 4.0_defReal-TWO*eps]
    u = [-ONE, ZERO, HALF]
    u = u/norm2(u)
    hs = surf % halfspace(r, u)
    if (.not.hs) then ! Perform small movment
      d = surf % distance(r, u)
      @assertTrue( abs(d) < 1.0E-6)
      @assertTrue(surf % halfspace(r + d*u, u))
    end if

  end subroutine testEdgeCases


end module box_test
