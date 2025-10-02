module wedge_test
  use numPrecision
  use universalVariables
  use dictionary_class,  only : dictionary
  use wedge_class,       only : wedge
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
  !! wedge Test case
  !!
  @TestCase(constructor=newTestCase)
    type, extends(ParameterizedTestCase) :: test_wedge
      integer(shortInt)                  :: axis
      integer(shortInt), dimension(2)    :: plane
      type(wedge)                        :: surf
    contains
      procedure :: tearDown
    end type test_wedge

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
         string = 'xWedge'
       case(Y_AXIS)
         string = 'yWedge'
       case(Z_AXIS)
         string = 'zWedge'
       case default
         string ="Unknown"
      end select

  end function toString

  !!
  !! Build new test_wedge test case
  !! Given integer direction X_AXIS, Y_AXIS or Z_AXIS
  !!
  function newTestCase(dir) result(tst)
    type(dirParam), intent(in)  :: dir
    type(test_wedge)            :: tst
    type(dictionary)            :: dict
    character(nameLen)          :: type
    real(defReal), dimension(3) :: origin

    ! Select type of wedge and axis
    select case(dir % dir)
      case(X_AXIS)
        tst % axis = X_AXIS
        tst % plane = [Y_AXIS, Z_AXIS]
        type = 'xWedge'

      case(Y_AXIS)
        tst % axis = Y_AXIS
        tst % plane = [X_AXIS, Z_AXIS]
        type = 'yWedge'

      case(Z_AXIS)
        tst % axis = Z_AXIS
        tst % plane = [X_AXIS, Y_AXIS]
        type = 'zWedge'

      case default
        print *, "Should not happen. Wrong direction in testCase constructor"

    end select

    origin = ZERO
    origin(tst % axis)     = 5.0_defReal
    origin(tst % plane(1)) = ONE

    ! Build surface
    call dict % init(7)
    call dict % store('id', 8)
    call dict % store('type', type)
    call dict % store('origin', origin)
    call dict % store('halfwidth', 4.0_defReal)
    call dict % store('height', 3.0_defReal)
    call dict % store('opening', 30.0_defReal)
    call dict % store('rotation', 60.0_defReal)
    call tst % surf % init(dict)

  end function newTestCase

  !!
  !! Deconstruct the test case
  !!
  subroutine tearDown(this)
    class(test_wedge), intent(inout) :: this

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
    class(test_wedge), intent(inout) :: this
    real(defReal), dimension(6)      :: aabb, ref
    character(nameLen)               :: name
    real(defReal), parameter         :: TOL = 1.0E-6_defReal

    ! Test ID
    @assertEqual(8, this % surf % id())

    ! Change ID
    call this % surf % setID(21)
    @assertEqual(21, this % surf % id())

    ! Name
    select case(this % axis)
      case(X_AXIS)
        name = 'xWedge'
      case(Y_AXIS)
        name = 'yWedge'
      case(Z_AXIS)
        name = 'zWedge'
    end select
    @assertEqual(name, this % surf % myType())

    ! Bounding Box
    ref(this % axis)      = ONE
    ref(this % axis + 3)  = 9.0_defReal
    ref(this % plane(1))  = ONE
    ref(this % plane(2))  = ZERO
    ref(this % plane(1) + 3) = 4.0_defReal
    ref(this % plane(2) + 3) = 6.0_defReal / sqrt(3.0_defReal)
    aabb = this % surf % boundingBox()

    @assertEqual(ref, aabb, TOL)

    ! Tolerance
    @assertEqual(this % surf % surfTol(), SURF_TOL * 4.0_defReal, TOL)

  end subroutine testMisc

  !!
  !! Test boundary conditions
  !!
! @Test(cases=[1,2,3])
!   subroutine testBC(this)
!     class(test_wedge), intent(inout) :: this
!     real(defReal), dimension(3)          :: r, u, r_pre, u_pre
!
!     ! Set boundary conditions
!     ! Should ignore extra entries
!     call this % surf % setBC([VACUUM_BC, REFLECTIVE_BC, REFLECTIVE_BC])
!
!     ! Apply BC
!     r = [ONE, ONE, ONE]
!     u = ZERO
!     ! Moving out at the surface in one planar direction
!     u(this % plane(1)) = ONE
!
!     r_pre = r
!     u_pre = u
!
!     ! Explicit
!     call this % surf % explicitBC(r, u)
!     @assertEqual(r_pre, r)
!     @assertEqual(u_pre, u)
!
!     ! Transform
!     call this % surf % transformBC(r, u)
!     @assertEqual(r_pre, r)
!     @assertEqual(u_pre, u)
!
!   end subroutine testBC

  !!
  !! Test Halfspaces membership
  !!
@Test(cases=[1,2,3])
  subroutine testHalfspace(this)
    class(test_wedge), intent(inout) :: this
    integer(shortInt)                :: a, p1, p2
    real(defReal), dimension(3)      :: r, u
    real(defReal)                    :: eps, tolerance

    ! Get surface tolerance
    tolerance = this % surf % surfTol()

    ! Set axis and plane axis indices
    a  = this % axis
    p1 = this % plane(1)
    p2 = this % plane(2)

    ! Obviously out
    r(a)  = 7.0_defReal
    r(p1) = 3.0_defReal
    r(p2) = 6.0_defReal
    u(a)  = ZERO
    u(p1) = ONE
    u(p2) = ZERO
    @assertTrue(this % surf % halfspace(r, u))

    ! Obviously in
    r(a)  = 7.0_defReal
    r(p1) = 3.0_defReal
    r(p2) = sqrt(3.0_defReal)
    @assertFalse(this % surf % halfspace(r, u))

    ! Choose point at the surface of the wedge moving in
    r(p1) = ONE
    u(a)  = TWO
    u = u / norm2(u)
    @assertFalse(this % surf % halfspace(r, u))

    ! Choose point at the surface of the wedge moving out
    u(p1) = -ONE
    u = u / norm2(u)
    @assertTrue(this % surf % halfspace(r, u))

    ! Out within tolerance going out
    eps = 0.01_defReal * tolerance
    @assertTrue(this % surf % halfspace(r + eps*u, u))

    ! Out within tolerance going in
    @assertFalse(this % surf % halfspace(r + eps*u, -u))

    ! Out outside tolerance
    eps = TWO * tolerance
    @assertTrue(this % surf % halfspace(r + eps*u, u))

    ! Well outside
    eps = HALF
    @assertTrue(this % surf % halfspace(r + eps*u, u))

    ! Well within
    eps = -HALF
    @assertFalse(this % surf % halfspace(r + eps*u, u))

    ! In within tolerance going in
    eps = -0.01_defReal * tolerance
    @assertFalse(this % surf % halfspace(r + eps*u, -u))

    ! Tangent particle should be outside
    u     = ZERO
    u(p2) = -ONE
    u = u /norm2(u)
    @assertTrue(this % surf % halfspace(r - eps*u, u))

    ! ! Choose point on one of the bases of the wedge, moving out
    r(a)       = 9.0_defReal
    r([p1,p2]) = TWO
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
    class(test_wedge), intent(inout) :: this
    integer(shortInt)                :: a, p1, p2
    real(defReal), dimension(3)      :: r, u
    real(defReal)                    :: ref, tolerance
    real(defReal), parameter :: TOL = 1.0E-7

    ! Get surface tolerance
    tolerance = this % surf % surfTol()

    ! Set axis and plane axis indices
    a = this % axis
    p1 = this % plane(1)
    p2 = this % plane(2)

    ! **Outside the wedge
    r(a)  = 4.0_defReal
    r(p1) = ZERO
    r(p2) = TWO

    ! Impact at an angle, direction of flight on an axis
    u(a)  = ONE / sqrt(TWO)
    u(p1) = ONE / sqrt(TWO)
    u(p2) = ZERO
    ref   = sqrt(TWO)
    @assertEqual(ref, this % surf % distance(r, u), TOL * ref)

    ! Impacting the bases
    r(a) = -ONE
    ref  = TWO * sqrt(TWO)
    @assertEqual(ref, this % surf % distance(r, u), TOL * ref)

    ! Parallel outside
    u(a)  = ONE
    u(p1) = ZERO
    u(p2) = ZERO
    u   = u/norm2(u)
    @assertEqual(INF, this % surf % distance(r, u))

    ! Going in the opposite direction
    u(a)  = ONE
    u(p1) = -ONE
    u(p2) = -ONE
    u   = u/norm2(u)
    @assertEqual(INF, this % surf % distance(r, u))

    ! Parallel to the bases
    u(a)  = ZERO
    u(p1) = 3.0_defReal
    u(p2) = 4.0_defReal
    u   = u/norm2(u)
    @assertEqual(INF, this % surf % distance(r, u))

    ! Parallel to the bases going in
    r(a)  = 4.0_defReal
    r(p1) = HALF * (7.0_defReal + sqrt(3.0_defReal))
    r(p2) = HALF * (5.0_defReal * sqrt(3.0_defReal) - ONE)
    u(a)  = ZERO
    u(p1) = -HALF
    u(p2) = -HALF * sqrt(3.0_defReal)
    u   = u/norm2(u)
    ref = TWO
    @assertEqual(ref, this % surf % distance(r, u))

    ! **Exactly at the surface
    r(p1) = ONE
    r(p2) = 3.0_defReal * HALF / sqrt(3.0_defReal)

    ! Particle going inside
    u(a)  = ZERO
    u(p1) = ONE
    u(p2) = ZERO
    u = u/norm2(u)
    ref = 1.5_defReal
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
    r(p1) = ONE - 0.1_defReal * tolerance
    u(a)  = ZERO
    u(p1) = ONE
    u(p2) = ZERO
    u = u/norm2(u)
    ref = 1.5_defReal + 0.1_defReal * tolerance
    @assertEqual(ref, this % surf % distance(r, u), TOL * ref)

    ! **Inside within surface tolerance
    r(p1) = ONE + 0.1_defReal * tolerance
    u(a)  = -ONE
    u(p1) = ZERO
    u(p2) = ZERO
    u = u/norm2(u)
    @assertEqual(INF, this % surf % distance(r, u))

    ! ** Inside the surface
    r(a)  = 7.0_defReal
    r(p1) = 3.0_defReal
    r(p2) = sqrt(3.0_defReal)

    ! Parallel to axis
    u(a)  = -ONE
    u(p1) = ZERO
    u(p2) = ZERO
    ref = 6.0_defReal
    @assertEqual(ref, this % surf % distance(r, u), TOL * ref)

    ! Parallel to bases
    u(a)  = ZERO
    u(p1) = -ONE
    u(p2) = ZERO
    u = u/norm2(u)
    ref = TWO
    @assertEqual(ref, this % surf % distance(r, u), TOL * ref)

    ! Hitting the lower base
    u(a)  = -3.0_defReal
    u(p1) = HALF
    u = u/norm2(u)
    ref = sqrt(37.0_defReal)
    @assertEqual(ref, this % surf % distance(r, u), TOL * ref)

    ! ** Outside the surface and on top
    r(a)  = 10.0_defReal
    r(p1) = 6.0_defReal
    r(p2) = ONE

    ! Parallel to axis
    u(a)  = -ONE
    u(p1) = ZERO
    u(p2) = ZERO
    @assertEqual(INF, this % surf % distance(r, u))

    ! At an angle
    u(p1) = -0.1_defReal
    u = u/norm2(u)
    @assertEqual(INF, this % surf % distance(r, u))

    ! At an angle parallel to wedge
    u(a)  = ZERO
    u(p1) = -ONE
    u(p2) = ONE
    u = u/norm2(u)
    @assertEqual(INF, this % surf % distance(r, u))

  end subroutine testDistance

end module wedge_test
