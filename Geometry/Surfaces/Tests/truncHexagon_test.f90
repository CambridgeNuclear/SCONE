module truncHexagon_test
  use numPrecision
  use universalVariables
  use dictionary_class,   only : dictionary
  use truncHexagon_class, only : truncHexagon
  use funit

  implicit none

  !!
  !! Test parameter wrapper around TWO INTEGERS (bit of boilerplate)
  !!
  @testParameter(constructor=newParam)
  type, extends (AbstractTestParameter) :: hexParam
     integer(shortInt) :: dir
     integer(shortInt) :: orient
  contains
     procedure :: toString
  end type hexParam

  !!
  !! Truncated Hexagon Test case
  !!
  @TestCase(constructor=newTestCase)
    type, extends(ParameterizedTestCase) :: test_truncHexagon
      integer(shortInt)               :: axis
      integer(shortInt), dimension(2) :: plane
      integer(shortInt)               :: orient
      type(truncHexagon)              :: surf
    contains
      procedure :: tearDown
    end type test_truncHexagon

contains

  !!
  !! Test parameter constructor
  !!
  function newParam(i) result(param)
     integer(shortInt), intent(in) :: i
     type (hexParam) :: param

     param % orient = mod(i,10)
     param % dir = i / 10

  end function newParam

  !!
  !! Print parameter to string for more verbose description
  !!
  function toString(this) result(string)
     class (hexParam), intent(in) :: this
     character(:), allocatable :: string

     select case(this % dir)
       case(X_AXIS)
         string = 'xTruncHexagon'
       case(Y_AXIS)
         string = 'yTruncHexagon'
       case(Z_AXIS)
         string = 'zTruncHexagon'
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
  !! Halfwidth            3
  !! Halfheight 5.0
  !! Orientation          1
  !! ID 75
  !!
  function newTestCase(param) result(tst)
    type(hexParam), intent(in)  :: param
    type(test_truncHexagon)     :: tst
    type(dictionary)            :: dict
    character(nameLen)          :: type
    real(defReal), dimension(3) :: origin
    real(defReal)               :: hw, hh
    integer(shortInt)           :: orient

    ! Select type of squareCylinder and axis
    orient = param % orient
    select case(param % dir)
      case(X_AXIS)
        tst % axis = X_AXIS
        if (orient == 1) then
          tst % plane = [Y_AXIS, Z_AXIS]
        else
          tst % plane = [Z_AXIS, Y_AXIS]
        end if
        type = 'xTruncHexagon'

      case(Y_AXIS)
        tst % axis = Y_AXIS
        if (orient == 1) then
          tst % plane = [X_AXIS, Z_AXIS]
        else
          tst % plane = [Z_AXIS, X_AXIS]
        end if
        type = 'yTruncHexagon'

      case(Z_AXIS)
        tst % axis = Z_AXIS
        if (orient == 1) then
          tst % plane = [X_AXIS, Y_AXIS]
        else
          tst % plane = [Y_AXIS, X_AXIS]
        end if
        type = 'zTruncHexagon'

      case default
        print *, "Should not happen. Wrong direction in testcase constructor"
        error stop

    end select

    ! Set origin, halfwidth and orientation
    origin(tst % plane(1)) = ONE
    origin(tst % plane(2)) = TWO
    origin(tst % axis) = TWO

    hw = 3.0_defReal
    hh = 5.0_defReal

    ! Build surface
    call dict % init(6)
    call dict % store('id', 75)
    call dict % store('type', type)
    call dict % store('origin', origin)
    call dict % store('halfwidth', hw)
    call dict % store('halfheight', hh)
    call dict % store('orientation', orient)
    call tst % surf % init(dict)

  end function newTestCase

  !!
  !! Deconstruct the test case
  !!
  subroutine tearDown(this)
    class(test_truncHexagon), intent(inout) :: this

    call this % surf % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Proper tests begin here
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    ! Points on the flat will have the coordinates:
    ! (+4/-1, [2-3sqrt(3)/4, 2+3sqrt(3)/4])
    !
    ! Points on the top right have the coordinates:
    ! x+sqrt(3)y=7+2sqrt(3) (for x in [1,4])
    !
    ! Points on the top left have the coordinates:
    ! -x+sqrt(3)y=5+2sqrt(3) (for x in [-2,1])
    !
    ! Points on the bottom left have the coordinates:
    ! -x-sqrt(3)y=5-2sqrt(3) (for x in [-2,1])
    !
    ! Points on the bottom right have the coordinates:
    ! x-sqrt(3)y=7-2sqrt(3) (for x in [1,4])
    !
    ! Span along y is 2/sqrt(3) * halfwidth

  !!
  !! Test Misc functionality
  !!
  !! Directions must be given as integers for pFUnit parser to work
  !!
@Test(cases = [11, 12, 21, 22, 31, 32])
  subroutine testMisc(this)
    class(test_truncHexagon), intent(inout) :: this
    real(defReal), dimension(6)             :: aabb, ref
    character(nameLen)                      :: name
    real(defReal), parameter                :: TOL = 1.0E-6_defReal

    ! Test ID
    @assertEqual(75, this % surf % id())

    ! Change ID
    call this % surf % setID(1)
    @assertEqual(1, this % surf % id())

    ! Name
    select case(this % axis)
      case(X_AXIS)
        name = 'xTruncHexagon'
      case(Y_AXIS)
        name = 'yTruncHexagon'
      case(Z_AXIS)
        name = 'zTruncHexagon'
    end select
    @assertEqual(name, this % surf % myType())

    ! Bounding Box
    ! For orientation 1, span along x is 2 * halfwidth.
    ! Span along y is 2/sqrt(3) * halfwidth
    ref = [ZERO, ZERO, ZERO, 4.0_defReal, 4.0_defReal, 4.0_defReal]
    ref(this % plane(1)) = -2.0_defReal
    ref(this % plane(2)) = TWO - TWO/sqrt(3.0_defReal) * 3.0_defReal
    ref(this % plane(1)+3) = 4.0_defReal
    ref(this % plane(2)+3) = TWO + TWO/sqrt(3.0_defReal) * 3.0_defReal
    ref(this % axis) = -3.0_defReal
    ref(this % axis+3) = 7.0_defReal

    aabb = this % surf % boundingBox()
    @assertEqual(ref, aabb, TOL)

  end subroutine testMisc

  !!
  !! Test boundary conditions
  !!
@Test(cases = [11, 12, 21, 22, 31, 32])
  subroutine testBC(this)
    class(test_truncHexagon), intent(inout) :: this
    integer(shortInt), dimension(6)         :: BC
    integer(shortInt)                       :: ax, p1, p2
    real(defReal), dimension(3)             :: r, u, r_ref, u_ref, n, periodicShift
    real(defReal), parameter                :: TOL = 1.0E-6

    ! Get axis and diffrent planar directions
    ax = this % axis
    p1 = this % plane(1)
    p2 = this % plane(2)

    ! Explicit BCs

    ! Vacuum test
    BC = VACUUM_BC
    call this % surf % setBC(BC)
    
    ! Travelling out the left edge
    r([ax, p1, p2]) = [ZERO, -TWO, 3.0_defReal]
    u([ax, p1, p2]) = [ZERO, -ONE, ZERO]
    r_ref = r
    u_ref = u
    call this % surf % explicitBC(r, u)
    @assertEqual(r_ref, r, TOL)
    @assertEqual(u_ref, u, TOL)
    
    ! Travelling out the right edge
    r([ax, p1, p2]) = [ZERO, 4.0_defReal, 3.0_defReal]
    u([ax, p1, p2]) = [ZERO, ONE, ZERO]
    r_ref = r
    u_ref = u
    call this % surf % explicitBC(r, u)
    @assertEqual(r_ref, r, TOL)
    @assertEqual(u_ref, u, TOL)

    ! Travelling out the top right-inclined edge
    r([ax, p1, p2]) = [ZERO, 2.0_defReal, 5.0_defReal/sqrt(3.0_defReal) + TWO]
    u([ax, p1, p2]) = [ZERO, ONE/TWO, sqrt(3.0_defReal)/TWO]
    r_ref = r
    u_ref = u
    call this % surf % explicitBC(r, u)
    @assertEqual(r_ref, r, TOL)
    @assertEqual(u_ref, u, TOL)

    ! Travelling out the bottom
    r([ax, p1, p2]) = [-3.0_defReal, 2.0_defReal, 2.0_defReal]
    u([ax, p1, p2]) = [-ONE, ZERO, ZERO]
    r_ref = r
    u_ref = u
    call this % surf % explicitBC(r, u)
    @assertEqual(r_ref, r, TOL)
    @assertEqual(u_ref, u, TOL)

    ! Axial reflection test
    BC(ax) = REFLECTIVE_BC
    BC(ax + 3) = REFLECTIVE_BC
    call this % surf % setBC(BC)
   
    r([ax, p1, p2]) = [-3.0_defReal, 2.0_defReal, 2.0_defReal]
    u([ax, p1, p2]) = [-ONE, ZERO, ZERO]
    r_ref = r
    u_ref([ax, p1, p2]) = [ONE, ZERO, ZERO]
    call this % surf % explicitBC(r, u)
    @assertEqual(r_ref, r, TOL)
    @assertEqual(u_ref, u, TOL)

    ! Reflection test
    !BC = REFLECTIVE_BC
    !call this % surf % setBC(BC)
    
    ! Travelling out the right edge
    !r([ax, p1, p2]) = [ZERO, 4.0_defReal, 3.0_defReal]
    !u([ax, p1, p2]) = [ZERO, ONE, ZERO]
    !r_ref = r
    !u_ref = u
    !u_ref(p1) = -ONE
    !call this % surf % explicitBC(r, u)
    !@assertEqual(r_ref, r, TOL)
    !@assertEqual(u_ref, u, TOL)

    ! Travelling out the bottom right edge
    ! x-sqrt(3)y=7-2sqrt(3) (for x in [1,4])
    ! Normal is i/2 - j*sqrt(3)/2
    !r([ax, p1, p2]) = [ZERO, 3.0_defReal, -4.0_defReal/sqrt(3.0_defReal)+TWO]
    !u([ax, p1, p2]) = [ZERO, sqrt(TWO)/TWO, -sqrt(TWO)/TWO]
    !r_ref = r
    !n = [ZERO, HALF, -sqrt(3.0_defReal)/2]
    !u_ref = u - TWO * dot_product(n, u) * n
    !call this % surf % explicitBC(r, u)
    !@assertEqual(r_ref, r, TOL)
    !@assertEqual(u_ref, u, TOL)

    ! Periodic test
    BC = PERIODIC_BC
    call this % surf % setBC(BC)

    ! Leave the bottom left
    ! -x-sqrt(3)y=5-2sqrt(3) (for x in [-2,1])
    ! Normal is -i/2 -j * sqrt(3)/2
    r([ax, p1, p2]) = [ZERO, ZERO, -5.0_defReal/sqrt(3.0_defReal) + TWO]
    u([ax, p1, p2]) = [ZERO, ZERO, ONE]
    periodicShift([ax, p1, p2]) = [ZERO, 3.0_defReal, sqrt(3.0_defReal) * 3.0_defReal]
    r_ref([ax, p1, p2]) = r([ax, p1, p2]) + periodicShift([ax, p1, p2])
    u_ref = u
    call this % surf % explicitBC(r, u)
    @assertEqual(r_ref, r, TOL)
    @assertEqual(u_ref, u, TOL)

    ! Corner
    ! Should arrive at the bottom left due to direction
    r([ax, p1, p2]) = [ZERO, ONE, TWO + TWO * sqrt(3.0_defReal)]
    u([ax, p1, p2]) = [ZERO, ONE, ZERO]
    u = u/norm2(u)

    r_ref([ax, p1, p2]) = [ZERO, -TWO, TWO - sqrt(3.0_defReal)]
    u_ref = u
    
    call this % surf % explicitBC(r, u)
    @assertEqual(r_ref, r, TOL)
    @assertEqual(u_ref, u, TOL)

    ! Top plane
    r([ax, p1, p2]) = [7.0_defReal, 2.0_defReal, 2.0_defReal]
    u([ax, p1, p2]) = [ONE, ZERO, ZERO]
    r_ref([ax, p1, p2]) = [-3.0_defReal, 2.0_defReal, 2.0_defReal]
    u_ref([ax, p1, p2]) = [ONE, ZERO, ZERO]
    call this % surf % explicitBC(r, u)
    @assertEqual(r_ref, r, TOL)
    @assertEqual(u_ref, u, TOL)

    ! Transform BCs

    ! Periodic
    ! Shift the particle 5.25 lattices right
    r([ax, p1, p2]) = [ZERO, 32.5_defReal, 1.9_defReal]
    u([ax, p1, p2]) = [ZERO, ONE, -ONE]
    u = u/norm2(u)

    r_ref([ax, p1, p2]) = [ZERO, 2.5_defReal, 1.9_defReal]
    u_ref = u

    call this % surf % transformBC(r, u)
    @assertEqual(r_ref, r, TOL)
    @assertEqual(u_ref, u, TOL)

    ! Try the corner
    ! Shift to the lattice above the point
    r([ax, p1, p2]) = [10.0_defReal, ONE, TWO + sqrt(3.0_defReal) * 6.0_defReal]
    u([ax, p1, p2]) = [ZERO, ONE, -ONE]
    u = u/norm2(u)

    r_ref([ax, p1, p2]) = [ZERO, ONE, TWO]
    u_ref([ax, p1, p2]) = [ZERO, ONE, -ONE]
    u_ref = u_ref/norm2(u_ref)

    call this % surf % transformBC(r, u)
    @assertEqual(r_ref, r, TOL)
    @assertEqual(u_ref, u, TOL)

    ! Reflective
    ! Axial only
    BC = VACUUM_BC
    BC(ax) = REFLECTIVE_BC
    BC(ax + 3) = REFLECTIVE_BC
    call this % surf % setBC(BC)

    r([ax, p1, p2]) = [-4.0_defReal, TWO, TWO]
    u([ax, p1, p2]) = [-ONE, ONE, ZERO]
    u = u/norm2(u)

    r_ref([ax, p1, p2]) = [-TWO, TWO, TWO]
    u_ref([ax, p1, p2]) = [ONE, ONE, ZERO]
    u_ref = u_ref/norm2(u_ref)

    call this % surf % transformBC(r, u)
    @assertEqual(r_ref, r, TOL)
    @assertEqual(u_ref, u, TOL)

    !BC = REFLECTIVE_BC
    !call this % surf % setBC(BC)

    ! Shift left by one lattice starting just off centre
    !r([ax, p1, p2]) = [ZERO, -5.5_defReal, 1.5_defReal]
    !u([ax, p1, p2]) = [ZERO, ONE, -ONE]
    !u = u/norm2(u)
    
    !r_ref([ax, p1, p2]) = [ZERO, 1.5_defReal, 1.5_defReal]
    !u_ref([ax, p1, p2]) = [ZERO, -ONE, -ONE]
    !u_ref = u_ref/norm2(u_ref)

    !call this % surf % transformBC(r, u)
    !@assertEqual(r_ref, r, TOL)
    !@assertEqual(u_ref, u, TOL)
    
    ! Shift left by two lattice starting just off centre
    !r([ax, p1, p2]) = [ZERO, -10.5_defReal, 1.5_defReal]
    !u([ax, p1, p2]) = [ZERO, ONE, -ONE]
    !u = u/norm2(u)
    
    !r_ref([ax, p1, p2]) = [ZERO, 1.5_defReal, 1.5_defReal]
    !u_ref([ax, p1, p2]) = [ZERO, ONE, -ONE]
    !u_ref = u_ref/norm2(u_ref)

    !call this % surf % transformBC(r, u)
    !@assertEqual(r_ref, r, TOL)
    !@assertEqual(u_ref, u, TOL)

    ! Shift from the centre out the bottom right by just over
    ! half a lattice
    ! x-sqrt(3)y=7-2sqrt(3) (for x in [1,4])
    ! Normal is i/2 - j*sqrt(3)/2
    !r([ax, p1, p2]) = [ZERO, ONE + TWO, TWO - TWO * sqrt(3.0_defReal)]
    !u([ax, p1, p2]) = [ZERO, ONE, -ONE]
    !u = u/norm2(u)
    
    !r_ref([ax, p1, p2]) = [ZERO, TWO, TWO - sqrt(3.0_defReal)]
    !u_ref([ax, p1, p2]) = [ZERO, ONE - sqrt(3.0_defReal), ONE - sqrt(3.0_defReal)]
    !u_ref = u_ref/norm2(u_ref)

    !call this % surf % transformBC(r, u)
    !@assertEqual(r_ref, r, TOL)
    !@assertEqual(u_ref, u, TOL)

  end subroutine testBC

  !!
  !! Test Halfspaces membership
  !!
@Test(cases = [11, 12, 21, 22, 31, 32])
  subroutine testHalfspace(this)
    class(test_truncHexagon), intent(inout) :: this
    integer(shortInt)                       :: ax, p1, p2
    real(defReal), dimension(3)             :: r, u, u2
    real(defReal)                           :: eps

    ! Get axis and diffrent planar directions
    ax = this % axis
    p1 = this % plane(1)
    p2 = this % plane(2)

    ! ** Well inside the hexagon
    r = ZERO
    u = ZERO
    u(p2) = ONE
    r([p1, p2]) = [1.8_defReal, ONE]
    @assertFalse(this % surf % halfspace(r, u))

    ! Different Quadrant
    r([p1, p2]) = [-0.2_defReal, 2.1_defReal]
    @assertFalse(this % surf % halfspace(r, u))

    ! ** Well outside the hexagon
    r([p1, p2]) = [-4.0_defReal, 2.0_defReal]
    @assertTrue(this % surf % halfspace(r, u))

    ! Diffrent Quadrant
    r([p1, p2]) = [5.0_defReal, -3.0_defReal]
    @assertTrue(this % surf % halfspace(r, u))

    ! Above the hexagon
    r([p1, p2]) = [TWO, TWO]
    r(ax) = 12.0_defReal
    @assertTrue(this % surf % halfspace(r, u))

    ! Below
    r(ax) = -10.0_defReal
    @assertTrue(this % surf % halfspace(r, u))

    ! Back in
    r(ax) = 2.0_defReal
    @assertFalse(this % surf % halfspace(r, u))

    ! ** Proximity of the surface
    r([p1, p2]) = [-TWO, TWO]
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
    u = ZERO
    u(p1) = ONE

    @assertTrue(this % surf % halfspace(r + eps*u, u2))
    @assertTrue(this % surf % halfspace(r - eps*u, u2))
    
    ! Crossing the top surface
    r([ax, p1, p2]) = [7.0_defReal, TWO, TWO]
    u = ZERO
    u(ax) = -ONE
    @assertFalse(this % surf % halfspace(r, u))
    @assertTrue(this % surf % halfspace(r, -u))

    ! Within surface surface tolerance
    eps = -HALF * SURF_TOL
    @assertTrue(this % surf % halfspace(r + eps*u, -u))

    ! A bit above Surface Tolerance
    eps = 1.00001_defReal * SURF_TOL
    @assertFalse(this % surf % halfspace(r + eps*u, -u))

  end subroutine testHalfspace

  !!
  !! Test distance calculation
  !!
@Test(cases = [11, 12, 21, 22, 31, 32])
  subroutine testDistance(this)
    class(test_truncHexagon), intent(inout) :: this
    integer(shortInt)                       :: ax, p1, p2
    real(defReal), dimension(3)             :: r, u
    real(defReal)                           :: ref
    real(defReal), parameter :: TOL = 1.0E-7

    ! Get axis and different planar directions
    ax = this % axis
    p1 = this % plane(1)
    p2 = this % plane(2)

    !** Outside the hexagon
    r = ZERO
    r(p1) = -3.0_defReal
    r(p2) = TWO

    ! Direct hit
    u([ax, p1, p2]) = [ONE, ONE, ZERO]
    u = u/norm2(u)
    ref = SQRT2
    @assertEqual(ref, this % surf % distance(r, u), ref * TOL)

    ! Moving away
    @assertEqual(INF, this % surf % distance(r, -u))

    ! Moving straight up
    u = ZERO
    u(ax) = ONE
    @assertEqual(INF, this % surf % distance(r, u))

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
    r([ax, p1, p2]) = [ZERO, ZERO, TWO - TWO * sqrt(3.0_defReal) - SURF_TOL]
    u([ax, p1, p2])  = [ZERO, ONE, ZERO]
    u = u/norm2(u)
    @assertEqual(INF, this % surf % distance(r, u))
    
    r([ax, p1, p2]) = [ZERO, ZERO, TWO - TWO * sqrt(3.0_defReal) + SURF_TOL]
    u([ax, p1, p2])  = [ZERO, ONE, ZERO]
    u = u/norm2(u)
    @assertEqual(ONE, this % surf % distance(r, u), TOL)

    ! Below, hitting the bottom
    r([ax, p1, p2]) = [-6.0_defReal, TWO, TWO]
    u = ZERO
    u(ax) = ONE
    @assertEqual(3.0_defReal, this % surf % distance(r, u), TOL)

    ! And missing the bottom
    u(ax) = -ONE
    @assertEqual(INF, this % surf % distance(r, u), TOL)

    ! Parallel
    ! Should hit the bottom right slanted surface
    r([ax, p1, p2]) = [ZERO, -TWO, TWO - HALF * sqrt(3.0_defReal)]
    u([ax, p1, p2]) = [ZERO, HALF * sqrt(3.0_defReal), -HALF]
    u = u/norm2(u)
    @assertEqual(5.0_defReal * sqrt(3.0_defReal) * HALF, this % surf % distance(r, u), TOL)

    ! ** At the surface
    r([ax, p1, p2]) = [ZERO, ZERO, TWO - 5.0_defReal/sqrt(3.0_defReal)]
    u([ax, p1, p2]) = [ZERO, HALF, HALF * sqrt(3.0_defReal)]
    ref = 6.0_defReal
    @assertEqual(ref, this % surf % distance(r, u), ref * TOL)

    ! Outside within Surface Tolerance
    r([ax, p1, p2]) = [ZERO, ZERO, TWO - 5.0_defReal/sqrt(3.0_defReal)]
    u([ax, p1, p2]) = [ZERO, HALF, HALF * sqrt(3.0_defReal)]
    r([ax, p1, p2]) = r([ax, p1, p2]) - HALF * SURF_TOL * u([ax, p1, p2])
    ref = 6.0_defReal + HALF * SURF_TOL
    @assertEqual(ref, this % surf % distance(r, u), ref * TOL)

    ! Inside within Surface Tolerance
    r([ax, p1, p2]) = [ZERO, ZERO, TWO - 5.0_defReal/sqrt(3.0_defReal)]
    u([ax, p1, p2]) = [ZERO, HALF, HALF * sqrt(3.0_defReal)]
    r([ax, p1, p2]) = r([ax, p1, p2]) - HALF * SURF_TOL * u([ax, p1, p2])
    ref = 6.0_defReal - HALF * SURF_TOL
    @assertEqual(ref, this % surf % distance(r, u), ref * TOL)

    ! ** Inside
    r([ax, p1, p2]) = [ZERO, -1.8_defReal, TWO - sqrt(3.0_defReal) / 4]
    u([ax, p1, p2]) = [ZERO, HALF, HALF * sqrt(3.0_defReal)]
    ref = 79.0_defReal / 20
    @assertEqual(ref, this % surf % distance(r, u), ref * TOL)

    ! Other direction
    ref = 0.4_defReal
    @assertEqual(ref, this % surf % distance(r, -u), ref * TOL)

    ! Hitting the top
    u([ax, p1, p2]) = [ONE, ZERO, ZERO]
    @assertEqual(7.0_defReal, this % surf % distance(r, u), ref * TOL)
    @assertEqual(3.0_defReal, this % surf % distance(r, -u), ref * TOL)

  end subroutine testDistance

  !!
  !! Test normal calculation
  !!
@Test(cases = [11, 12, 21, 22, 31, 32])
  subroutine testNormal(this)
    class(test_truncHexagon), intent(inout) :: this
    integer(shortInt)                       :: ax, p1, p2
    real(defReal), dimension(3)             :: r, u, n
    real(defReal), parameter :: TOL = 1.0E-7
    
    ! Get axis and different planar directions
    ax = this % axis
    p1 = this % plane(1)
    p2 = this % plane(2)

    u = [8.0_defReal, -42.0_defReal, 13.1_defReal]

    ! Test on different hex faces
    r(p1) = -TWO
    r(p2) = TWO
    r(ax) = TWO

    n = this % surf % normal(r, u)
    @assertEqual(-ONE, n(p1), TOL)
    @assertEqual(ZERO, n(p2), TOL)
    @assertEqual(ZERO, n(ax), TOL)
    
    r(p1) = 4.0_defReal
    r(p2) = TWO
    r(ax) = TWO

    n = this % surf % normal(r, u)
    @assertEqual(ONE, n(p1), TOL)
    @assertEqual(ZERO, n(p2), TOL)
    @assertEqual(ZERO, n(ax), TOL)

    r(p1) = TWO
    r(p2) = TWO + 5.0_defReal/sqrt(3.0_defReal)
    
    n = this % surf % normal(r, u)
    @assertEqual(HALF, n(p1), TOL)
    @assertEqual(sqrt(3.0_defReal)*HALF, n(p2), TOL)
    @assertEqual(ZERO, n(ax), TOL)

    r(p1) = ZERO
    r(p2) = TWO - 5.0_defReal/sqrt(3.0_defReal)
    
    n = this % surf % normal(r, u)
    @assertEqual(-HALF, n(p1), TOL)
    @assertEqual(-sqrt(3.0_defReal)*HALF, n(p2), TOL)
    @assertEqual(ZERO, n(ax), TOL)

    ! Put on corners
    r(p1) = ONE
    r(p2) = TWO + TWO * sqrt(3.0_defReal)
    
    n = this % surf % normal(r, u)
    @assertEqual(ZERO, n(p1), TOL)
    @assertEqual(ONE, n(p2), TOL)
    @assertEqual(ZERO, n(ax), TOL)
    
    r(p1) = ONE
    r(p2) = TWO - TWO * sqrt(3.0_defReal)
    
    n = this % surf % normal(r, u)
    @assertEqual(ZERO, n(p1), TOL)
    @assertEqual(-ONE, n(p2), TOL)
    @assertEqual(ZERO, n(ax), TOL)

    r(p1) = -TWO
    r(p2) = TWO - sqrt(3.0_defReal)
    
    n = this % surf % normal(r, u)
    @assertEqual(-sqrt(3.0_defReal) * HALF, n(p1), TOL)
    @assertEqual(-HALF, n(p2), TOL)
    @assertEqual(ZERO, n(ax), TOL)
    
    ! Put on the top and bottom planes
    r(p1) = ONE
    r(p2) = TWO
    r(ax) = 7.0_defReal
    
    n = this % surf % normal(r, u)
    @assertEqual(ZERO, n(p1), TOL)
    @assertEqual(ZERO, n(p2), TOL)
    @assertEqual(ONE, n(ax), TOL)
    
    r(ax) = -3.0_defReal
    n = this % surf % normal(r, u)
    @assertEqual(ZERO, n(p1), TOL)
    @assertEqual(ZERO, n(p2), TOL)
    @assertEqual(-ONE, n(ax), TOL)

    ! And way below the bottom plane
    r(ax) = -30.0_defReal
    n = this % surf % normal(r, u)
    @assertEqual(ZERO, n(p1), TOL)
    @assertEqual(ZERO, n(p2), TOL)
    @assertEqual(-ONE, n(ax), TOL)
    
  end subroutine testNormal

end module truncHexagon_test
