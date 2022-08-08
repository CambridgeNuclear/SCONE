module cone_test
  use numPrecision
  use universalVariables
  use dictionary_class,  only : dictionary
  use cone_class,        only : cone
  use pfUnit_mod

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
    type, extends(ParameterizedTestCase) :: test_cone
      integer(shortInt)                  :: axis
      integer(shortInt), dimension(2)    :: plane
      type(cone)                         :: surf
    contains
      procedure :: tearDown
    end type test_cone

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
         string = 'xCone'
       case(Y_AXIS)
         string = 'yCone'
       case(Z_AXIS)
         string = 'zCone'
       case default
         string ="Unknown"
      end select
  end function toString

  !!
  !! Build new test_cone test case
  !! Given integer direction X_AXIS, Y_AXIS or Z_AXIS
  !!
  !! Vertex 1.0, 1.0, 1.0
  !! Tangent 1.0
  !! ID 52
  !!
  function newTestCase(dir) result(tst)
    type(dirParam), intent(in) :: dir
    type(test_cone)            :: tst
    type(dictionary)           :: dict
    character(nameLen)         :: type

    ! Select type of cone and axis
    select case(dir % dir)
      case(X_AXIS)
        tst % axis = X_AXIS
        tst % plane = [Y_AXIS, Z_AXIS]
        type = 'xCone'

      case(Y_AXIS)
        tst % axis = Y_AXIS
        tst % plane = [X_AXIS, Z_AXIS]
        type = 'yCone'

      case(Z_AXIS)
        tst % axis = Z_AXIS
        tst % plane = [X_AXIS, Y_AXIS]
        type = 'zCone'

      case default
        print *, "Should not happen. Wrong direction in testCase constructor"

    end select

    ! Build surface
    call dict % init(5)
    call dict % store('id', 52)
    call dict % store('type', type)
    call dict % store('vertex', [ONE, ONE, ONE])
    call dict % store('tangent', ONE)
    call dict % store('orientation', -1)
    call tst % surf % init(dict)
  end function newTestCase

  !!
  !! Deconstruct the test case
  !!
  subroutine tearDown(this)
    class(test_cone), intent(inout) :: this

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
    class(test_cone), intent(inout) :: this
    real(defReal), dimension(6)     :: aabb, ref
    character(nameLen)              :: name
    real(defReal), parameter        :: TOL = 1.0E-6_defReal

    ! Test ID
    @assertEqual(52, this % surf % id())

    ! Change ID
    call this % surf % setID(1)
    @assertEqual(1, this % surf % id())

    ! Name
    select case(this % axis)
      case(X_AXIS)
        name = 'xCone'
      case(Y_AXIS)
        name = 'yCone'
      case(Z_AXIS)
        name = 'zCone'
    end select
    @assertEqual(name, this % surf % myType())

    ! Bounding Box
    ref = [-INF, -INF, -INF, INF, INF, INF]
    ref( 3 + this % axis) = ONE
    aabb = this % surf % boundingBox()
    @assertEqual(ref, aabb, TOL)

  end subroutine testMisc

  !!
  !! Test boundary conditions
  !!
@Test(cases=[1,2,3])
  subroutine testBC(this)
    class(test_cone), intent(inout) :: this
    real(defReal), dimension(3)     :: r, u, r_pre, u_pre

    ! Set Boundary Contidions
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
    class(test_cone), intent(inout) :: this
    integer(shortInt)               :: a, p1, p2
    real(defReal), dimension(3)     :: r, u
    real(defReal)                   :: eps

    ! Set axis and plane axis indices
    a = this % axis
    p1 = this % plane(1)
    p2 = this % plane(2)

    ! Choose point at the surface moving in
    ! Direction is at 45deg to the plane
    r = ONE
    r(a) = ZERO
    r(p1) = TWO
    u(p2) = ZERO
    u(p1) = -ONE
    u(a)  = -ONE
    u = u / norm2(u)

    ! At the surface
    @assertFalse(this % surf % halfspace(r, u))

    ! Out within SURF_TOL
    eps = -HALF * SURF_TOL
    @assertFalse(this % surf % halfspace(r + eps*u, u))

    ! Out outside SURF_TOL
    eps = -TWO * SURF_TOL
    @assertTrue(this % surf % halfspace(r + eps*u, u))

    ! Well Outside
    eps = -TWO
    @assertTrue(this % surf % halfspace(r + eps*u, u))

    ! Well withn
    eps = TWO
    @assertFalse(this % surf % halfspace(r + eps*u, u))

    ! Tangent particle should be outside
    u = ZERO
    u(p2) = -ONE
    u = u /norm2(u)
    @assertTrue(this % surf % halfspace(r, u))

  end subroutine testHalfspace

  !!
  !! Test distance calculation
  !!
@Test(cases=[1, 2, 3])
  subroutine testDistance(this)
    class(test_cone), intent(inout) :: this
    integer(shortInt)               :: a, p1, p2
    real(defReal), dimension(3)     :: r, u
    real(defReal)                   :: ref
    real(defReal), parameter :: TOL = 1.0E-7

    ! Set axis and plane axis indices
    a = this % axis
    p1 = this % plane(1)
    p2 = this % plane(2)

    ! **Outside the cone
    r(a) = ZERO
    r(p1) = 3.0_defReal
    r(p2) = ONE

    ! Impact at an angle, direction of flight on an axis
    u(a) = ZERO
    u(p1) = -ONE
    u(p2) = ZERO
    u = u/norm2(u)
    ref = ONE
    @assertEqual(ref, this % surf % distance(r, u), TOL * ref)

    ! Almost Parallel
    u(a) = -ONE
    u(p1) = 1.0_defReal + 1.0E+20_defReal
    u(p2) = ZERO
    u = u/norm2(u)
    @assertEqual(INF, this % surf % distance(r, u))

    ! **Exactly at the surface
    r(p1) = TWO

    ! Particle going inside
    u(a) = ZERO
    u(p1) = -ONE
    u(p2) = ZERO
    u = u/norm2(u)
    ref = TWO
    @assertEqual(ref, this % surf % distance(r, u), TOL * ref)

    ! Particle going outside
    @assertEqual(INF, this % surf % distance(r, -u))

    ! Tangent particle
    u(a) = ZERO
    u(p1) = ZERO
    u(p2) = ONE
    u = u/norm2(u)
    @assertEqual(INF, this % surf % distance(r, u))

    ! **Outside within surface tolerance
    r(p1) = TWO + HALF * SURF_TOL
    u(a) = ZERO
    u(p1) = -ONE
    u(p2) = ZERO
    u = u/norm2(u)
    ref = TWO + HALF * SURF_TOL
    @assertEqual(ref, this % surf % distance(r, u), TOL * ref)

    ! ** Inside the surface
    r = ONE
    r(a) = ZERO

    ! Parallel to second plane direction
    u(a) = ZERO
    u(p1) = ZERO
    u(p2) = ONE
    u = u/norm2(u)
    ref = ONE
    @assertEqual(ref, this % surf % distance(r, u), TOL * ref)

    ! Parallel to axis
    u(a) = ONE
    u(p1) = ZERO
    u(p2) = ZERO
    u = u/norm2(u)
    ref = ONE
    @assertEqual(ref, this % surf % distance(r, u), TOL * ref)

    u(a) = -ONE
    u = u/norm2(u)
    @assertEqual(INF, this % surf % distance(r, u), TOL * ref)

    ! At an angle
    u(a) = ONE
    u(p1) = ONE
    u(p2) = ZERO
    u = u/norm2(u)
    ref = HALF * SQRT2

    @assertEqual(ref, this % surf % distance(r, u), TOL * ref)

  end subroutine testDistance

end module cone_test
