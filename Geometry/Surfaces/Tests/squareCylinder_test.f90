module squareCylinder_test
  use numPrecision
  use universalVariables
  use dictionary_class,     only : dictionary
  use squareCylinder_class, only : squareCylinder
  use funit

  implicit none

  !!
  !! Test parameter wrapper around AN INTEGER (bit of boilerplate)
  !!
  !!
  @testParameter(constructor=newParam)
  type, extends (AbstractTestParameter) :: dirParam
     integer(shortInt) :: dir
  contains
     procedure :: toString
  end type dirParam

  !!
  !! Cylinder Test case
  !!
  @TestCase(constructor=newTestCase)
    type, extends(ParameterizedTestCase) :: test_squareCylinder
      integer(shortInt)               :: axis
      integer(shortInt), dimension(2) :: plane
      type(squareCylinder)            :: surf
    contains
      procedure :: tearDown
    end type test_squareCylinder

contains

  !!
  !! Test parameter constructor
  !!
  function newParam(i) result(param)
     integer(shortInt), intent(in) :: i
     type (dirParam) :: param

     param % dir = i

  end function newParam

  !!
  !! Print parameter to string for more verbose description
  !!
  function toString(this) result(string)
     class (dirParam), intent(in) :: this
     character(:), allocatable :: string

     select case(this % dir)
       case(X_AXIS)
         string = 'xSquareCylinder'
       case(Y_AXIS)
         string = 'ySquareCylinder'
       case(Z_AXIS)
         string = 'zSquareCylinder'
       case default
         string ="Unknown"
      end select
  end function toString

  !!
  !! Build new test_squareCylinder test case
  !! Given integer direction X_AXIS, Y_AXIS or Z_AXIS
  !!
  !!            axis p1   p2
  !! Origin     2.0, 1.0, 2.0
  !! Halfwidths      2.0, 3.0
  !! ID 75
  !!
  function newTestCase(dir) result(tst)
    type(dirParam), intent(in) :: dir
    type(test_squareCylinder)  :: tst
    type(dictionary)           :: dict
    character(nameLen)         :: type
    real(defReal), dimension(3) :: origin, hw

    ! Select type of squareCylinder and axis
    select case(dir % dir)
      case(X_AXIS)
        tst % axis = X_AXIS
        tst % plane = [Y_AXIS, Z_AXIS]
        type = 'xSquareCylinder'

      case(Y_AXIS)
        tst % axis = Y_AXIS
        tst % plane = [X_AXIS, Z_AXIS]
        type = 'ySquareCylinder'

      case(Z_AXIS)
        tst % axis = Z_AXIS
        tst % plane = [X_AXIS, Y_AXIS]
        type = 'zSquareCylinder'

      case default
        print *, "Should not happen. Wrong direction in testcase constructor"

    end select

    ! Set origin & halfwidth
    origin = TWO
    origin(tst % plane(1)) = ONE
    origin(tst % plane(2)) = TWO

    hw = ZERO
    hw(tst % plane(1)) = TWO
    hw(tst % plane(2)) = 3.0_defReal

    ! Build surface
    call dict % init(4)
    call dict % store('id', 75)
    call dict % store('type', type)
    call dict % store('origin', origin)
    call dict % store('halfwidth', hw)
    call tst % surf % init(dict)

  end function newTestCase

  !!
  !! Deconstruct the test case
  !!
  subroutine tearDown(this)
    class(test_squareCylinder), intent(inout) :: this

    call this % surf % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Proper tests begin here
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test Misc functionality
  !!
  !! Directions must be given as integers for pFUnit parser to work
  !!
@Test(cases = [1, 2, 3])
  subroutine testMisc(this)
    class(test_squareCylinder), intent(inout) :: this
    real(defReal), dimension(6) :: aabb, ref
    character(nameLen)          :: name
    real(defReal), parameter    :: TOL = 1.0E-6_defReal

    ! Test ID
    @assertEqual(75, this % surf % id())

    ! Change ID
    call this % surf % setID(1)
    @assertEqual(1, this % surf % id())

    ! Name
    select case(this % axis)
      case(X_AXIS)
        name = 'xSquareCylinder'
      case(Y_AXIS)
        name = 'ySquareCylinder'
      case(Z_AXIS)
        name = 'zSquareCylinder'
    end select
    @assertEqual(name, this % surf % myType())

    ! Bounding Box
    ref = [ZERO, ZERO, ZERO, 4.0_defReal, 4.0_defReal, 4.0_defReal]
    ref(this % plane) = [-1.0_defReal, -1.0_defReal]
    ref(this % plane+3) = [3.0_defReal, 5.0_defReal]
    ref(this % axis) = -INF
    ref(this % axis+3) = INF

    aabb = this % surf % boundingBox()
    @assertEqual(ref, aabb, TOL)

  end subroutine testMisc

  !!
  !! Test boundary conditions
  !!
@Test(cases=[1,2,3])
  subroutine testBC(this)
    class(test_squareCylinder), intent(inout) :: this
    integer(shortInt), dimension(6)  :: BC
    integer(shortInt)                :: ax, p1, p2
    real(defReal), dimension(3)      :: r, u, r_ref, u_ref
    real(defReal), parameter :: TOL = 1.0E-6

    ! Get axis and diffrent planar directions
    ax = this % axis
    p1 = this % plane(1)
    p2 = this % plane(2)

    ! Set BC
    BC = VACUUM_BC
    BC(p1*2) = REFLECTIVE_BC
    BC(p1*2-1) = VACUUM_BC
    BC(p2*2) = PERIODIC_BC
    BC(p2*2-1) = PERIODIC_BC

    call this % surf % setBC(BC)

    ! Explicit BC

    ! Vacuum face
    r([ax, p1, p2]) = [ZERO, -1.0_defReal, 3.0_defReal]
    u([ax, p1, p2]) = [ZERO, -ONE, ZERO]
    r_ref = r
    u_ref = u
    call this % surf % explicitBC(r, u)
    @assertEqual(r_ref, r, TOL)
    @assertEqual(u_ref, u, TOL)

    ! Reflection face
    r([ax, p1, p2]) = [ZERO, 3.0_defReal, 3.0_defReal]
    u([ax, p1, p2]) = [ZERO, ONE, ZERO]
    r_ref = r
    u_ref = u
    u_ref(p1) = -ONE
    call this % surf % explicitBC(r, u)
    @assertEqual(r_ref, r, TOL)
    @assertEqual(u_ref, u, TOL)

    ! Periodic face
    r([ax, p1, p2]) = [ZERO, 2.0_defReal, 5.0_defReal]
    u([ax, p1, p2]) = [ZERO, ZERO, ONE]
    r_ref([ax, p1, p2]) = [ZERO, 2.0_defReal, -1.0_defReal]
    u_ref = u
    call this % surf % explicitBC(r, u)
    @assertEqual(r_ref, r, TOL)
    @assertEqual(u_ref, u, TOL)

    ! Corner
    r([ax, p1, p2]) = [ZERO, 3.0_defReal, -1.0_defReal]
    u([ax, p1, p2]) = [ZERO, ONE, -ONE]
    u = u/norm2(u)
    r_ref([ax, p1, p2]) = [ZERO, 3.0_defReal, 5.0_defReal]
    u_ref([ax, p1, p2]) = [ZERO, -ONE, -ONE]
    u_ref = u_ref/norm2(u_ref)
    call this % surf % explicitBC(r, u)
    @assertEqual(r_ref, r, TOL)
    @assertEqual(u_ref, u, TOL)

    ! Transform BC
    r([ax, p1, p2]) = [ZERO, 20.0_defReal, -13.0_defReal]
    u([ax, p1, p2]) = [ZERO, ONE, -ONE]
    u = u/norm2(u)

    r_ref([ax, p1, p2]) = [ZERO, -14.0_defReal, 5.0_defReal]
    u_ref([ax, p1, p2]) = [ZERO, -ONE, -ONE]
    u_ref = u_ref/norm2(u_ref)

    call this % surf % transformBC(r, u)
    @assertEqual(r_ref, r, TOL)
    @assertEqual(u_ref, u, TOL)

    ! Try the corner
    r([ax, p1, p2]) = [ZERO, 3.0_defReal, -1.0_defReal]
    u([ax, p1, p2]) = [ZERO, ONE, -ONE]
    u = u/norm2(u)

    r_ref([ax, p1, p2]) = [ZERO, 3.0_defReal, 5.0_defReal]
    u_ref([ax, p1, p2]) = [ZERO, -ONE, -ONE]
    u_ref = u_ref/norm2(u_ref)

    call this % surf % transformBC(r, u)
    @assertEqual(r_ref, r, TOL)
    @assertEqual(u_ref, u, TOL)

  end subroutine testBC

  !!
  !! Test Halfspaces membership
  !!
@Test(cases=[1,2,3])
  subroutine testHalfspace(this)
    class(test_squareCylinder), intent(inout) :: this
    integer(shortInt)                         :: ax, p1, p2
    real(defReal), dimension(3)               :: r, u, u2
    real(defReal)                             :: eps

    ! Get axis and diffrent planar directions
    ax = this % axis
    p1 = this % plane(1)
    p2 = this % plane(2)

    ! ** Well inside the square cylinder
    r = ZERO
    u = ZERO
    u(p2) = ONE
    r([p1, p2]) = [2.0_defReal, 0.0_defReal]
    @assertFalse(this % surf % halfspace(r, u))

    ! Different Quadrant
    r([p1, p2]) = [-0.5_defReal, 2.0_defReal]
    @assertFalse(this % surf % halfspace(r, u))

    ! ** Well outside the square cylinder
    r([p1, p2]) = [-1.5_defReal, 2.0_defReal]
    @assertTrue(this % surf % halfspace(r, u))

    ! Diffrent Quadrant
    r([p1, p2]) = [0.5_defReal, 5.2_defReal]
    @assertTrue(this % surf % halfspace(r, u))

    ! ** Proximity of the surface
    r([p1, p2]) = [-1.0_defReal, 3.0_defReal]
    u = ZERO
    u(p1) = -ONE
    @assertTrue(this % surf % halfspace(r, u))
    @assertFalse(this % surf % halfspace(r, -u))

    ! Inside within Surface Tolerance
    eps = -HALF * SURF_TOL
    @assertTrue(this % surf % halfspace(r + eps*u, u))

    ! A bit above Surface Tolerance
    eps = -1.0001_defReal * SURF_TOL
    @assertFalse(this % surf % halfspace(r + eps*u, u))

    ! Parallel to the surface
    u2 = ZERO
    u2(p2) = ONE
    eps = HALF * SURF_TOL

    @assertTrue(this % surf % halfspace(r + eps*u, u2))
    @assertFalse(this % surf % halfspace(r - eps*u, u2))

    ! At the surface in diffrent quadrant



  end subroutine testHalfspace

  !!
  !! Test distance calculation
  !!
@Test(cases=[1, 2, 3])
  subroutine testDistance(this)
    class(test_squareCylinder), intent(inout) :: this
    integer(shortInt)                   :: ax, p1, p2
    real(defReal), dimension(3)         :: r, u
    real(defReal)                       :: ref
    real(defReal), parameter :: TOL = 1.0E-7

    ! Get axis and diffrent planar directions
    ax = this % axis
    p1 = this % plane(1)
    p2 = this % plane(2)

    !** Outside the square cylinder
    r = ZERO
    r(p1) = -2.0_defReal
    r(p2) = ZERO

    ! Direct hit
    u([ax, p1, p2]) = [ONE, ONE, ZERO]
    u = u/norm2(u)
    ref = SQRT2
    @assertEqual(ref, this % surf % distance(r, u), ref * TOL)

    ! Moving away
    @assertEqual(INF, this % surf % distance(r, -u))

    ! Oblique Hit
    u([ax, p1, p2]) = [ONE, ONE, ONE]
    u = u/norm2(u)
    ref = sqrt(3.0_defReal)
    @assertEqual(ref, this % surf % distance(r, u), ref * TOL)

    ! Miss
    u([ax, p1, p2]) = [ONE, ONE, -TWO]
    u = u/norm2(u)
    @assertEqual(INF, this % surf % distance(r, u))

    ! Corner skim
    ! Use dirty values
    r([ax, p1, p2]) = [ZERO, -1.3_defReal, 0.3_defReal]
    u([ax, p1, p2])  = [ZERO, 0.3_defReal, -1.3_defReal]
    u = u/norm2(u)
    @assertEqual(INF, this % surf % distance(r, u))

    ! Parallel
    r([ax, p1, p2]) = [ZERO, -1.3_defReal, 0.3_defReal]
    u([ax, p1, p2])  = [ZERO, ZERO, ONE]
    u = u/norm2(u)
    @assertEqual(INF, this % surf % distance(r, u))

    ! ** At the surface
    r([ax, p1, p2]) = [ZERO, -1.0_defReal, 0.5_defReal]
    u([ax, p1, p2]) = [ZERO, ONE, ZERO]
    ref = 4.0_defReal
    @assertEqual(ref, this % surf % distance(r, u), ref * TOL)

    ! Outside within Surface Tolerance
    r([ax, p1, p2]) = [ZERO, -1.0_defReal - HALF * SURF_TOL, 0.5_defReal]
    u([ax, p1, p2]) = [ZERO, ONE, ZERO]
    ref = 4.0_defReal + HALF * SURF_TOL
    @assertEqual(ref, this % surf % distance(r, u), ref * TOL)

    ! Inside within Surface Tolerance
    r([ax, p1, p2]) = [ZERO, -1.0_defReal + HALF * SURF_TOL, 0.5_defReal]
    u([ax, p1, p2]) = [ZERO, ONE, ZERO]
    ref = 4.0_defReal - HALF * SURF_TOL
    @assertEqual(ref, this % surf % distance(r, u), ref * TOL)

    ! ** Inside
    r([ax, p1, p2]) = [9.0_defReal, 1.15_defReal, 1.0_defReal]
    u([ax, p1, p2]) = [ZERO, ZERO, ONE]
    ref = 4.0_defReal
    @assertEqual(ref, this % surf % distance(r, u), ref * TOL)

    ! Other direction
    ref = 2.0_defReal
    @assertEqual(ref, this % surf % distance(r, -u), ref * TOL)

  end subroutine testDistance

  !!
  !! Test Edge Cases
  !!
  !! Test unlikley cases to make sure that halfspace + distance
  !! procedures allow particle to escape
  !!
@Test(cases=[1, 2, 3])
  subroutine testEdgeCases(this)
    class(test_squareCylinder), intent(inout) :: this
    real(defReal), dimension(3) :: r, u, u2
    integer(shortInt)           :: ax, p1, p2
    real(defReal)               :: eps, d
    logical(defBool)            :: hs

    ! Get axis and diffrent planar directions
    ax = this % axis
    p1 = this % plane(1)
    p2 = this % plane(2)

    ! ** Corner
    ! * Particle is almost at the corner
    !   Either it is outside or can escape with a short movment in next step
    eps =  5.0_defReal * epsilon(eps)
    r([ax, p1, p2]) = [ZERO, -1.0_defReal+eps, -1.0_defReal+eps]
    u([ax, p1, p2]) = [ZERO, TWO, -ONE]
    u = u/norm2(u)
    hs = this % surf % halfspace(r, u)
    if (.not.hs) then ! Perform small movment
      d = this % surf % distance(r, u)
      @assertTrue( abs(d) < 1.0E-6)
      @assertTrue(this % surf % halfspace(r + d*u, u))
    end if

    ! Try in other direction
    u([ax, p1, p2]) = [ZERO, -ONE, TWO]
    u = u/norm2(u)
    hs = this % surf % halfspace(r, u)
    if (.not.hs) then ! Perform small movment
      d = this % surf % distance(r, u)
      @assertTrue( abs(d) < 1.0E-6)
      @assertTrue(this % surf % halfspace(r + d*u, u))
    end if

    ! Try asymertic corner
    ! Point is not exactly at the diagonal
    r([ax, p1, p2]) = [ZERO, -1.0_defReal+eps, -1.0_defReal+eps*TWO]
    u([ax, p1, p2]) = [ZERO, TWO, -ONE]
    u = u/norm2(u)
    hs = this % surf % halfspace(r, u)
    if (.not.hs) then ! Perform small movment
      d = this % surf % distance(r, u)
      @assertTrue( abs(d) < 1.0E-6)
      @assertTrue(this % surf % halfspace(r + d*u, u))
    end if

    ! Try in other direction
    r([ax, p1, p2]) = [ZERO, -1.0_defReal+eps*TWO, -1.0_defReal+eps]
    u([ax, p1, p2]) = [ZERO, -ONE, TWO]
    u = u/norm2(u)
    hs = this % surf % halfspace(r, u)
    if (.not.hs) then ! Perform small movment
      d = this % surf % distance(r, u)
      @assertTrue( abs(d) < 1.0E-6)
      @assertTrue(this % surf % halfspace(r + d*u, u))
    end if

  end subroutine testEdgeCases

  !!
  !! Test encountered problems
  !!
  !! Contains test related to bugs found at some point
  !! TODO: Move some of this tests to main test procedures
  !!
@Test(cases=[1])
  subroutine test_problems(this)
    class(test_squareCylinder), intent(inout) :: this ! Ignore this
    type(squareCylinder) :: surf
    type(dictionary)     :: dict
    real(defReal), dimension(3) :: r, u

    call dict % init(5)
    call dict % store('type','zSquareCylinder')
    call dict % store('id', 7)
    call dict % store('origin', [ZERO, ZERO, ZERO])
    call dict % store('halfwidth', [8.00_defReal, 1.26_defReal, 0.0_defReal])
    call surf % init(dict)

    r = [-7.63_defReal, 1.26_defReal, ZERO]
    u = [ZERO, -ONE, ZERO]
    @assertFalse(surf % halfspace(r, u))

  end subroutine test_problems



end module squareCylinder_test
