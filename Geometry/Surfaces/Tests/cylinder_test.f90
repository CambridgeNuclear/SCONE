module cylinder_test
  use numPrecision
  use universalVariables
  use dictionary_class,  only : dictionary
  use cylinder_class,    only : cylinder
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
    type, extends(ParameterizedTestCase) :: test_cylinder
      integer(shortInt)               :: axis
      integer(shortInt), dimension(2) :: plane
      type(cylinder)                  :: surf
    contains
      procedure :: tearDown
    end type test_cylinder

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
         string = 'xCylinder'
       case(Y_AXIS)
         string = 'yCylinder'
       case(Z_AXIS)
         string = 'zCylinder'
       case default
         string ="Unknown"
      end select
  end function toString

  !!
  !! Build new test_cylinder test case
  !! Given integer direction X_AXIS, Y_AXIS or Z_AXIS
  !!
  !! Origin 2.0, 2.0, 2.0
  !! Radius 2.0
  !! ID 75
  !!
  function newTestCase(dir) result(tst)
    type(dirParam), intent(in) :: dir
    type(test_cylinder)        :: tst
    type(dictionary)      :: dict
    character(nameLen)    :: type

    ! Select type of cylinder and axis
    select case(dir % dir)
      case(X_AXIS)
        tst % axis = X_AXIS
        tst % plane = [Y_AXIS, Z_AXIS]
        type = 'xCylinder'

      case(Y_AXIS)
        tst % axis = Y_AXIS
        tst % plane = [X_AXIS, Z_AXIS]
        type = 'yCylinder'

      case(Z_AXIS)
        tst % axis = Z_AXIS
        tst % plane = [X_AXIS, Y_AXIS]
        type = 'zCylinder'

      case default
        print *, "Should not happen. Wrong direction in testcase constructor"

    end select

    ! Build surface
    call dict % init(4)
    call dict % store('id', 75)
    call dict % store('type', type)
    call dict % store('origin', [TWO, TWO, TWO])
    call dict % store('radius', TWO)
    call tst % surf % init(dict)
  end function newTestCase

  !!
  !! Deconstruct the test case
  !!
  subroutine tearDown(this)
    class(test_cylinder), intent(inout) :: this

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
    class(test_cylinder), intent(inout) :: this
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
        name = 'xCylinder'
      case(Y_AXIS)
        name = 'yCylinder'
      case(Z_AXIS)
        name = 'zCylinder'
    end select
    @assertEqual(name, this % surf % myType())

    ! Bounding Box
    ref = [ZERO, ZERO, ZERO, 4.0_defReal, 4.0_defReal, 4.0_defReal]
    ref(this % axis) = -INF
    ref( 3 + this % axis) = INF
    aabb = this % surf % boundingBox()
    @assertEqual(ref, aabb, TOL)

  end subroutine testMisc

  !!
  !! Test boundary conditions
  !!
@Test(cases=[1,2,3])
  subroutine testBC(this)
    class(test_cylinder), intent(inout) :: this
    real(defReal), dimension(3) :: r, u, r_pre, u_pre

    ! Set Boundary Contidions
    ! Should ignore extra entries
    call this % surf % setBC([VACUUM_BC, REFLECTIVE_BC, REFLECTIVE_BC])

    ! Apply BC
    r = [TWO, TWO, TWO]
    u = ZERO
    ! Moving out at the surface in one planar direction
    r(this % plane(1)) = 4.0_defReal
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
    class(test_cylinder), intent(inout) :: this
    integer(shortInt)                   :: a, p1, p2
    real(defReal), dimension(3)    :: r, u
    real(defReal)                  :: eps

    ! Set axis and plane axis indices
    a = this % axis
    p1 = this % plane(1)
    p2 = this % plane(2)

    ! Choose point at the surface moving in
    ! DIrection is at 45deg to the plane
    r = TWO
    r(p1) = ZERO
    u(p2) = ZERO
    u(p1) = ONE
    u(a)  = ONE
    u = u / norm2(u)

    ! At the surface
    @assertFalse(this % surf % halfspace(r, u))

    ! Out within SURF_TOL
    eps = -SURF_TOL
    @assertFalse(this % surf % halfspace(r + eps*u, u))

    ! Out outside SURF_TOL
    eps = -SQRT2 * 1.00001_defReal * SURF_TOL
    @assertTrue(this % surf % halfspace(r + eps*u, u))

    ! Well Outside
    eps = -TWO
    @assertTrue(this % surf % halfspace(r + eps*u, u))

    ! Well withn
    eps = TWO
    @assertFalse(this % surf % halfspace(r + eps*u, u))

    ! Tangent particle should be outside
    u(p2) =  ONE
    u(p1) = ZERO
    u(a)  = ONE
    u = u /norm2(u)
    @assertTrue(this % surf % halfspace(r, u))

  end subroutine testHalfspace

  !!
  !! Test distance calculation
  !!
@Test(cases=[1, 2, 3])
  subroutine testDistance(this)
    class(test_cylinder), intent(inout) :: this
    integer(shortInt)                   :: a, p1, p2
    real(defReal), dimension(3)         :: r, u
    real(defReal)                       :: ref
    real(defReal), parameter :: TOL = 1.0E-7

    ! Set axis and plane axis indices
    a = this % axis
    p1 = this % plane(1)
    p2 = this % plane(2)

    ! **Outside the cylinder
    r = TWO
    r(p1) = -ONE

    ! Perpendicular (in plane) impact
    u(a) = ONE
    u(p1) = ONE
    u(p2) = ZERO
    u = u/norm2(u)
    ref = SQRT2
    @assertEqual(ref, this % surf % distance(r, u), TOL * ref)

    ! Oblique impact
    ! at 30deg
    u(a)  = ONE
    u(p1) = sqrt(3.0_defReal) / TWO
    u(p2) = HALF
    u = u/norm2(u)
    ref = 1.275200556_defReal * SQRT2
    @assertEqual(ref, this % surf % distance(r, u), TOL * ref)

    ! Almost Parallel
    u(a) = 1.0E+20_defReal
    u(p1) = 1.0_defReal
    u(p2) = ZERO
    u = u/norm2(u)
    @assertEqual(INF, this % surf % distance(r, u))

    ! **Exactly at the surface
    r = TWO
    r(p1) = ZERO

    ! Particle going inside
    u(a) = ONE
    u(p1) = ONE
    u(p2) = ZERO
    u = u/norm2(u)
    ref = 4.0_defReal * SQRT2
    @assertEqual(ref, this % surf % distance(r, u), TOL * ref)

    ! Particle going outside
    @assertEqual(INF, this % surf % distance(r, -u))

    ! Tangent particle
    u(a) = ONE
    u(p1) = ZERO
    u(p2) = ONE
    u = u/norm2(u)
    @assertEqual(INF, this % surf % distance(r, -u))

    ! **Outside within surface tolerance
    r = TWO
    r(p1) = ZERO - HALF * SURF_TOL
    u(a) = ONE
    u(p1) = ONE
    u(p2) = ZERO
    u = u/norm2(u)
    ref = (4.0_defReal + HALF * SURF_TOL) * SQRT2
    @assertEqual(ref, this % surf % distance(r, u), TOL * ref)

    ! ** Inside the surface
    r = TWO
    r(p1) = ONE

    ! +ve direction
    u(a) = ONE
    u(p1) = ONE
    u(p2) = ZERO
    u = u/norm2(u)
    ref = 3.0_defReal * SQRT2
    @assertEqual(ref, this % surf % distance(r, u), TOL * ref)

    ! -ve direction
    ref = 1.0_defReal * SQRT2
    @assertEqual(ref, this % surf % distance(r, -u), TOL * ref)

  end subroutine testDistance

end module cylinder_test
