module hexagon_test
  use numPrecision
  use universalVariables
  use dictionary_class,  only : dictionary
  use hexagon_class,     only : hexagon, flatType, pointType
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
    type, extends(ParameterizedTestCase) :: test_hexagon
      integer(shortInt)                  :: axis
      integer(shortInt)                  :: orient
      integer(shortInt), dimension(2)    :: plane
      type(hexagon)                      :: surf
    contains
      procedure :: tearDown
    end type test_hexagon

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
         string = 'xHexagon'
       case(Y_AXIS)
         string = 'yHexagon'
       case(Z_AXIS)
         string = 'zHexagon'
       case default
         string ="Unknown"
      end select
  end function toString

  !!
  !! Build new test_hexagon test case
  !! Given integer direction X_AXIS, Y_AXIS or Z_AXIS
  !!
  !!            axis p1   p2
  !! Origin     2.0, 1.0, 2.0
  !! Halfwidths      2.0, 3.0
  !! ID 75
  !!
  function newTestCase(dir) result(tst)
    type(dirParam), intent(in) :: dir
    type(test_hexagon)         :: tst
    type(dictionary)           :: dict
    character(nameLen)         :: type
    real(defReal), dimension(3) :: origin, hw

    ! Select type of hexagon and axis
    select case(dir % dir)
      case(X_AXIS)
        tst % axis = X_AXIS
        tst % plane = [Y_AXIS, Z_AXIS]
        type = 'xHexagon'

      case(Y_AXIS)
        tst % axis = Y_AXIS
        tst % plane = [X_AXIS, Z_AXIS]
        type = 'yHexagon'

      case(Z_AXIS)
        tst % axis = Z_AXIS
        tst % plane = [X_AXIS, Y_AXIS]
        type = 'zHexagon'

      case default
        print *, "Should not happen. Wrong direction in testcase constructor"
        error stop

    end select

    ! Set origin & halfwidth
    origin = TWO
    origin(tst % plane(1)) = ONE
    origin(tst % plane(2)) = TWO

    hw = TWO

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
    class(test_hexagon), intent(inout) :: this

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
    class(test_hexagon), intent(inout) :: this
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
        name = 'xHexagon'
      case(Y_AXIS)
        name = 'yHexagon'
      case(Z_AXIS)
        name = 'zHexagon'
    end select
    @assertEqual(name, this % surf % myType())

    ! Bounding Box

  end subroutine testMisc

  !!
  !! Test boundary conditions
  !!
@Test(cases=[1,2,3])
  subroutine testBC(this)
    class(test_hexagon), intent(inout) :: this
    integer(shortInt), dimension(6)    :: BC
    integer(shortInt)                  :: ax, p1, p2
    real(defReal), dimension(3)        :: r, u, r_ref, u_ref
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

  end subroutine testBC

  !!
  !! Test Halfspaces membership
  !!
@Test(cases=[1,2,3])
  subroutine testHalfspace(this)
    class(test_hexagon), intent(inout) :: this
    integer(shortInt)                  :: ax, p1, p2
    real(defReal), dimension(3)        :: r, u, u2
    real(defReal)                      :: eps

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


  end subroutine testHalfspace

  !!
  !! Test distance calculation
  !!
@Test(cases=[1, 2, 3])
  subroutine testDistance(this)
    class(test_hexagon), intent(inout) :: this
    integer(shortInt)                  :: ax, p1, p2
    real(defReal), dimension(3)        :: r, u
    real(defReal)                      :: ref
    real(defReal), parameter :: TOL = 1.0E-7

    ! Get axis and different planar directions
    ax = this % axis
    p1 = this % plane(1)
    p2 = this % plane(2)

    !** Outside the hexagon
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

  end subroutine testDistance

  !!
  !! Test Edge Cases
  !!
  !! Test unlikley cases to make sure that halfspace + distance
  !! procedures allow particle to escape
  !!
@Test(cases=[1, 2, 3])
  subroutine testEdgeCases(this)
    class(test_hexagon), intent(inout) :: this
    real(defReal), dimension(3) :: r, u
    integer(shortInt)           :: ax, p1, p2
    real(defReal)               :: eps, d
    logical(defBool)            :: hs


  end subroutine testEdgeCases

  !!
  !! Test encountered problems
  !!
  !! Contains test related to bugs found at some point
  !! TODO: Move some of this tests to main test procedures
  !!
@Test(cases=[1])
  subroutine test_problems(this)
    class(test_hexagon), intent(inout) :: this ! Ignore this
    type(squareCylinder) :: surf
    type(dictionary)     :: dict
    real(defReal), dimension(3) :: r, u


  end subroutine test_problems



end module hexagon_test
