module aPlane_test
  use numPrecision
  use universalVariables
  use dictionary_class,  only : dictionary
  use aPlane_class,      only : aPlane
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
    type, extends(ParameterizedTestCase) :: test_aPlane
      integer(shortInt)             :: axis
      type(aPlane)                  :: surf
    contains
      procedure :: tearDown
    end type test_aPlane

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
         string = 'xPlane'
       case(Y_AXIS)
         string = 'yPlane'
       case(Z_AXIS)
         string = 'zPlane'
       case default
         string ="Unknown"
      end select
  end function toString

  !!
  !! Build new test_aPlane test case
  !! Given integer direction X_AXIS, Y_AXIS or Z_AXIS
  !!
  !! ID 75
  !! x0/y0/z0 4.3;
  !!
  function newTestCase(dir) result(tst)
    type(dirParam), intent(in) :: dir
    type(test_aPlane)        :: tst
    type(dictionary)      :: dict

    ! Start dictionary
    ! Build surface
    call dict % init(4)
    call dict % store('id', 75)

    ! Select type of cylinder and axis
    select case(dir % dir)
      case(X_AXIS)
        tst % axis = X_AXIS
        call dict % store('type','xPlane')
        call dict % store('x0',4.3_defReal)

      case(Y_AXIS)
        tst % axis = Y_AXIS
        call dict % store('type','yPlane')
        call dict % store('y0',4.3_defReal)

      case(Z_AXIS)
        tst % axis = Z_AXIS
        call dict % store('type','zPlane')
        call dict % store('z0',4.3_defReal)

      case default
        print *, "Should not happen. Wrong direction in testcase constructor"
        error stop

    end select

    ! Build surface
    call tst % surf % init(dict)

  end function newTestCase

  !!
  !! Deconstruct the test case
  !!
  subroutine tearDown(this)
    class(test_aPlane), intent(inout) :: this

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
    class(test_aPlane), intent(inout) :: this
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
        name = 'xPlane'
      case(Y_AXIS)
        name = 'yPlane'
      case(Z_AXIS)
        name = 'zPlane'
    end select
    @assertEqual(name, this % surf % myType())

    ! Bounding Box
    ref = [-INF, -INF, -INF, INF, INF, INF]
    ref(this % axis) = 4.3_defReal
    ref( 3 + this % axis) = 4.3_defReal
    aabb = this % surf % boundingBox()
    @assertEqual(ref, aabb, TOL)

  end subroutine testMisc

  !!
  !! Test boundary conditions
  !!
@Test(cases=[1,2,3])
  subroutine testBC(this)
    class(test_aPlane), intent(inout) :: this
    real(defReal), dimension(3) :: r, u, r_pre, u_pre

    ! Set Boundary Contidions
    ! Should ignore extra entries
    call this % surf % setBC([VACUUM_BC, REFLECTIVE_BC, REFLECTIVE_BC])

    ! Apply BC
    r = [TWO, TWO, TWO]
    u = ZERO

    ! Put atthe surface
    r(this % axis) = 4.3_defReal
    u(this % axis) = ONE

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
    class(test_aPlane), intent(inout) :: this
    real(defReal), dimension(3)       :: r, u, u2
    real(defReal)                     :: eps

    r = ZERO
    u = ZERO

    ! At the surface
    r(this % axis) = 4.3_defReal
    u(this % axis) = ONE
    @assertTrue(this % surf % halfspace(r, u))

    ! In within SURF_TOL
    eps = -HALF * SURF_TOL
    @assertTrue(this % surf % halfspace(r + eps*u, u))

    ! In outside SURF_TOL
    eps = -1.0001_defReal * SURF_TOL
    @assertFalse(this % surf % halfspace(r + eps*u, u))

    ! Well Outside
    eps = TWO
    @assertTrue(this % surf % halfspace(r + eps*u, u))

    ! Well withn
    eps = -TWO
    @assertFalse(this % surf % halfspace(r + eps*u, u))

    ! Tangent particle should use position
    eps = HALF * SURF_TOL
    u2 = cshift(u, 1) ! Point to an orthogonal direction

    @assertTrue(this % surf % halfspace(r + eps*u, u2))
    @assertFalse(this % surf % halfspace(r - eps*u, u2))

  end subroutine testHalfspace

  !!
  !! Test distance calculation
  !!
@Test(cases=[1, 2, 3])
  subroutine testDistance(this)
    class(test_aPlane), intent(inout) :: this
    real(defReal), dimension(3)         :: r, u, u2
    real(defReal)                       :: ref
    real(defReal), parameter :: SQRT3 = sqrt(3.0_defReal)
    real(defReal), parameter :: TOL = 1.0E-7

    ! ** Inside
    r = TWO
    u = ZERO
    r(this % axis) = ZERO
    u = ONE   ! 45deg to each axis
    u = u/norm2(u)

    ref = 4.3_defReal * SQRT3
    @assertEqual(ref, this % surf % distance(r, u), TOL * ref)

    ! ** Exactly at surface
    r(this % axis) = 4.3_defReal
    @assertEqual(INF, this % surf % distance(r, u))

    ! ** Within Surface tolerance
    r(this % axis) = 4.3_defReal - HALF * SURF_TOL
    @assertEqual(INF, this % surf % distance(r, u))

    ! ** Well outside
    ! +ve direction
    r(this % axis) = 5.0_defReal
    @assertEqual(INF, this % surf % distance(r, u))

    ! -ve direction
    ref = 0.7_defReal * SQRT3
    @assertEqual(ref, this % surf % distance(r, -u), TOL * ref)

    ! ** Parallel to plane
    r = TWO
    r(this % axis) = ZERO
    u2 = ZERO
    u2(this % axis) = ONE
    u2 = cshift(u2, 1) ! Point to an orthogonal direction
    @assertEqual(INF, this % surf % distance(r, u2))

  end subroutine testDistance
  
  !!
  !! Test producing the normal vector
  !!
@Test(cases=[1, 2, 3])
  subroutine testNormal(this)
    class(test_aPlane), intent(inout) :: this
    real(defReal), dimension(3)       :: n
    real(defReal), dimension(3)       :: r, u

    r = [99.92_defReal, -6.0_defReal, 4.0_defReal]
    u = [100, 200, 400]

    n = this % surf % normal(r, u)

    @assertEqual(ONE, n(this % axis))

  end subroutine testNormal

end module aPlane_test
