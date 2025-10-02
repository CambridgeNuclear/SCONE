module truncCone_test
  use numPrecision
  use universalVariables
  use dictionary_class,  only : dictionary
  use truncCone_class,   only : truncCone
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
  !! Cone Test case
  !!
  @TestCase(constructor=newTestCase)
    type, extends(ParameterizedTestCase) :: test_truncCone
      integer(shortInt)                  :: axis
      integer(shortInt), dimension(2)    :: plane
      type(truncCone)                    :: surf
    contains
      procedure :: tearDown
    end type test_truncCone

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
         string = 'xTruncCone'
       case(Y_AXIS)
         string = 'yTruncCone'
       case(Z_AXIS)
         string = 'zTruncCone'
       case default
         string ="Unknown"
      end select

  end function toString

  !!
  !! Build new test_truncCone test case
  !! Given integer direction X_AXIS, Y_AXIS or Z_AXIS
  !!
  !! Vertex 1.0, 1.0, 1.0
  !! Tangent 1.0
  !! ID 52
  !!
  function newTestCase(dir) result(tst)
    type(dirParam), intent(in) :: dir
    type(test_truncCone)       :: tst
    type(dictionary)           :: dict
    character(nameLen)         :: type

    ! Select type of cone and axis
    select case(dir % dir)
      case(X_AXIS)
        tst % axis = X_AXIS
        tst % plane = [Y_AXIS, Z_AXIS]
        type = 'xTruncCone'

      case(Y_AXIS)
        tst % axis = Y_AXIS
        tst % plane = [X_AXIS, Z_AXIS]
        type = 'yTruncCone'

      case(Z_AXIS)
        tst % axis = Z_AXIS
        tst % plane = [X_AXIS, Y_AXIS]
        type = 'zTruncCone'

      case default
        print *, "Should not happen. Wrong direction in testCase constructor"

    end select

    ! Build surface
    call dict % init(5)
    call dict % store('id', 52)
    call dict % store('type', type)
    call dict % store('vertex', [ONE, ONE, ONE])
    call dict % store('angle', 45.0_defReal)
    call dict % store('hMin', ONE)
    call dict % store('hMax', 10.0_defReal)
    call tst % surf % init(dict)

  end function newTestCase

  !!
  !! Deconstruct the test case
  !!
  subroutine tearDown(this)
    class(test_truncCone), intent(inout) :: this

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
    class(test_truncCone), intent(inout) :: this
    real(defReal), dimension(6)          :: aabb, ref
    character(nameLen)                   :: name
    real(defReal), parameter             :: TOL = 1.0E-6_defReal

    ! Test ID
    @assertEqual(52, this % surf % id())

    ! Change ID
    call this % surf % setID(1)
    @assertEqual(1, this % surf % id())

    ! Name
    select case(this % axis)
      case(X_AXIS)
        name = 'xTruncCone'
      case(Y_AXIS)
        name = 'yTruncCone'
      case(Z_AXIS)
        name = 'zTruncCone'
    end select
    @assertEqual(name, this % surf % myType())

    ! Bounding Box
    ref(this % axis)      = TWO
    ref(this % axis + 3)  = 11.0_defReal
    ref(this % plane)     = -9.0_defReal
    ref(this % plane + 3) = 11.0_defReal
    aabb = this % surf % boundingBox()
    @assertEqual(ref, aabb, TOL)

    ! Tolerance
    @assertEqual(this % surf % surfTol(), SURF_TOL*5.50_defReal, TOL)

  end subroutine testMisc

  !!
  !! Test boundary conditions
  !!
@Test(cases=[1,2,3])
  subroutine testBC(this)
    class(test_truncCone), intent(inout) :: this
    real(defReal), dimension(3)          :: r, u, r_pre, u_pre

    ! Set boundary conditions
    ! Should ignore extra entries
    call this % surf % setBC([VACUUM_BC, REFLECTIVE_BC, REFLECTIVE_BC])

    ! Apply BC
    r = [ONE, ONE, ONE]
    u = ZERO
    ! Moving out at the surface in one planar direction
    u(this % plane(1)) = ONE

    r_pre = r
    u_pre = u

    ! Explicit
    call this % surf % explicitBC(r, u)
    @assertEqual(r_pre, r)
    @assertEqual(u_pre, u)

    ! Transform
    call this % surf % transformBC(r, u)
    @assertEqual(r_pre, r)
    @assertEqual(u_pre, u)

  end subroutine testBC

  !!
  !! Test Halfspaces membership
  !!
@Test(cases=[1,2,3])
  subroutine testHalfspace(this)
    class(test_truncCone), intent(inout) :: this
    integer(shortInt)                    :: a, p1, p2
    real(defReal), dimension(3)          :: r, u
    real(defReal)                        :: eps, tolerance

    ! Get surface tolerance
    tolerance = this % surf % surfTol()

    ! Set axis and plane axis indices
    a  = this % axis
    p1 = this % plane(1)
    p2 = this % plane(2)

    ! Choose point at the surface of the cone moving in
    ! Direction is at 45deg to the plane
    r         = ONE
    r([a,p1]) = 5.0_defReal
    u     = ZERO
    u(a)  = ONE
    u(p1) = -ONE
    u = u / norm2(u)

    ! At the surface
    @assertFalse(this % surf % halfspace(r, u))

    ! Out within tolerance
    ! NOTE that the tolerance is very sensitive: if using eps = HALF * TOL, the
    ! particle is outside tolerance
    eps = -0.01_defReal * tolerance
    @assertFalse(this % surf % halfspace(r + eps*u, u))

    ! Out outside tolerance
    eps = -TWO * tolerance
    @assertTrue(this % surf % halfspace(r + eps*u, u))

    ! Well outside
    eps = -HALF
    @assertTrue(this % surf % halfspace(r + eps*u, u))

    ! Well within
    eps = HALF
    @assertFalse(this % surf % halfspace(r + eps*u, u))

    ! Tangent particle should be outside
    u     = ZERO
    u(p2) = -ONE
    u = u /norm2(u)
    @assertTrue(this % surf % halfspace(r, u))

    ! Choose point on one of the bases of the cone, moving out
    r(a)       = 11.0_defReal
    r([p1,p2]) = 4.0_defReal
    u(a)       = ONE
    u([p1,p2]) = -ONE
    u = u / norm2(u)

    ! At the surface
    @assertTrue(this % surf % halfspace(r, u))

    ! In within tolerance
    eps = -HALF * tolerance
    @assertTrue(this % surf % halfspace(r + eps*u, u))

    ! In outside tolerance
    eps = -TWO * tolerance
    @assertFalse(this % surf % halfspace(r + eps*u, u))

    ! Well inside
    eps = -HALF
    @assertFalse(this % surf % halfspace(r + eps*u, u))

    ! Well outside
    eps = HALF
    @assertTrue(this % surf % halfspace(r + eps*u, u))

    ! Tangent particle should be outside
    u     = ZERO
    u(p2) = -ONE
    u = u /norm2(u)
    @assertTrue(this % surf % halfspace(r, u))

  end subroutine testHalfspace

  !!
  !! Test distance calculation
  !!
@Test(cases=[1, 2, 3])
  subroutine testDistance(this)
    class(test_truncCone), intent(inout) :: this
    integer(shortInt)                    :: a, p1, p2
    real(defReal), dimension(3)          :: r, u
    real(defReal)                        :: ref, tolerance
    real(defReal), parameter :: TOL = 1.0E-7

    ! Get surface tolerance
    tolerance = this % surf % surfTol()

    ! Set axis and plane axis indices
    a = this % axis
    p1 = this % plane(1)
    p2 = this % plane(2)

    ! **Outside the cone
    r(a)  = ZERO
    r(p1) = 4.0_defReal
    r(p2) = ONE

    ! Impact at an angle, direction of flight on an axis
    u(a)  = ONE
    u(p1) = ZERO
    u(p2) = ZERO
    ref   = 4.0_defReal
    @assertEqual(ref, this % surf % distance(r, u), TOL * ref)

    ! Impacting the bases
    u(a) = TWO
    u(p1) = -3.0_defReal
    u(p2) = ZERO
    u   = u/norm2(u)
    ref = sqrt(13.0_defReal)
    @assertEqual(ref, this % surf % distance(r, u), TOL * ref)

    ! Parallel outside
    u(a)  = ONE
    u(p1) = ONE
    u(p2) = ZERO
    u   = u/norm2(u)
    @assertEqual(INF, this % surf % distance(r, u))

    ! Going in the opposite direction
    u(a)  = -ONE
    u(p1) = ZERO
    u(p2) = ZERO
    u   = u/norm2(u)
    @assertEqual(INF, this % surf % distance(r, u))

    ! Parallel to the bases
    u(a)  = ZERO
    u(p1) = ONE
    u(p2) = ZERO
    u   = u/norm2(u)
    @assertEqual(INF, this % surf % distance(r, u))

    ! **Exactly at the surface
    r(a) = 4.0_defReal

    ! Particle going inside
    u(a)  = ZERO
    u(p1) = -ONE
    u(p2) = ZERO
    u = u/norm2(u)
    ref = 6.0_defReal
    @assertEqual(ref, this % surf % distance(r, u), TOL * ref)

    ! Particle going outside
    @assertEqual(INF, this % surf % distance(r, -u))

    ! Tangent particle
    u(a)  = ZERO
    u(p1) = ZERO
    u(p2) = ONE
    u = u/norm2(u)
    @assertEqual(INF, this % surf % distance(r, u))

    ! **Outside within surface tolerance
    r(p1) = 4.0_defReal + 0.1_defReal * tolerance
    u(a)  = ZERO
    u(p1) = -ONE
    u(p2) = ZERO
    u = u/norm2(u)
    ref = 6.0_defReal + HALF * tolerance
    @assertEqual(ref, this % surf % distance(r, u), TOL * ref)

    ! **Inside within surface tolerance
    r(a)  = 4.0_defReal + 0.1_defReal * tolerance
    r(p1) = 4.0_defReal
    u(a)  = -ONE
    u(p1) = ZERO
    u(p2) = ZERO
    u = u/norm2(u)
    @assertEqual(INF, this % surf % distance(r, u))

    ! ** Inside the surface
    r(a)  = 8.0_defReal
    r(p1) = 3.0_defReal
    r(p2) = ONE

    ! Parallel to axis
    u(a)  = ONE
    u(p1) = ZERO
    u(p2) = ZERO
    u = u/norm2(u)
    ref = 3.0_defReal
    @assertEqual(ref, this % surf % distance(r, u), TOL * ref)

    ! Parallel to bases
    u(a)  = ZERO
    u(p1) = -ONE
    u(p2) = ZERO
    u = u/norm2(u)
    ref = 9.0_defReal
    @assertEqual(ref, this % surf % distance(r, u), TOL * ref)

    ! Hitting at an angle
    u(a)  = -ONE
    u(p1) = ZERO
    u = u/norm2(u)
    ref = 5.0_defReal
    @assertEqual(ref, this % surf % distance(r, u), TOL * ref)

    ! Hitting at an angle, parallel to cone
    u(a)  = ONE
    u(p1) = -ONE
    u = u/norm2(u)
    ref = sqrt(TWO) * 3.0_defReal
    @assertEqual(ref, this % surf % distance(r, u), TOL * ref)

    ! ** Inside, hitting the lower base
    r(p1) = HALF
    u(a)  = -ONE
    u(p1) = ZERO
    ref = 6.0_defReal
    @assertEqual(ref, this % surf % distance(r, u), TOL * ref)

    ! ** Outside the surface and on top
    r(a)  = 10.0_defReal
    r(p1) = 13.0_defReal
    r(p2) = ONE

    ! Parallel to axis
    u(a)  = -ONE
    u(p1) = ZERO
    u(p2) = ZERO
    u = u/norm2(u)
    @assertEqual(INF, this % surf % distance(r, u))

    ! Other direction
    u(a)  = ONE
    u(p1) = ZERO
    u(p2) = ZERO
    u = u/norm2(u)
    @assertEqual(INF, this % surf % distance(r, u))

    ! At an angle
    u(a)  = ONE
    u(p1) = -ONE - HALF
    u(p2) = ZERO
    u = u/norm2(u)
    @assertEqual(INF, this % surf % distance(r, u))

    ! At an angle parallel to cone
    u(a)  = ONE
    u(p1) = -ONE
    u(p2) = ZERO
    u = u/norm2(u)
    @assertEqual(INF, this % surf % distance(r, u))

    ! ** To the right of the cone
    r(a)  = 16.0_defReal
    r(p1) = 4.0_defReal
    r(p2) = ONE

    ! Going in at an angle
    u(a)  = -5.0_defReal
    u(p1) = -2.0_defReal
    u(p2) = ZERO
    u = u/norm2(u)
    ref = sqrt(29.0_defReal)
    @assertEqual(ref, this % surf % distance(r, u), TOL * ref)

  end subroutine testDistance

end module truncCone_test
