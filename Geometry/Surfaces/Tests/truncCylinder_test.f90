module truncCylinder_test
  use numPrecision
  use universalVariables
  use dictionary_class,     only : dictionary
  use truncCylinder_class,  only : truncCylinder
  use funit

  implicit none

  !!
  !! Test parameter wrapper around AN INTEGER (bit of boilerplate)
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
    type, extends(ParameterizedTestCase) :: test_truncCylinder
      integer(shortInt)               :: axis
      integer(shortInt), dimension(2) :: plane
      type(truncCylinder)            :: surf
    contains
      procedure :: tearDown
    end type test_truncCylinder

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
         string = 'xTruncCylinder'
       case(Y_AXIS)
         string = 'yTruncCylinder'
       case(Z_AXIS)
         string = 'zTruncCylinder'
       case default
         string ="Unknown"
      end select
  end function toString

  !!
  !! Build new test_truncCylinder test case
  !! Given integer direction X_AXIS, Y_AXIS or Z_AXIS
  !!
  !!            axis p1   p2
  !! Origin     2.0, 1.0, 2.0
  !! Halfwidth  1.5
  !! Radius     2.0
  !! ID 75
  !!
  function newTestCase(dir) result(tst)
    type(dirParam), intent(in)  :: dir
    type(test_truncCylinder)    :: tst
    type(dictionary)            :: dict
    character(nameLen)          :: type
    real(defReal), dimension(3) :: origin
    real(defReal)               :: hw, r

    ! Select type of truncCylinder and axis
    select case(dir % dir)
      case(X_AXIS)
        tst % axis = X_AXIS
        tst % plane = [Y_AXIS, Z_AXIS]
        type = 'xTruncCylinder'

      case(Y_AXIS)
        tst % axis = Y_AXIS
        tst % plane = [X_AXIS, Z_AXIS]
        type = 'yTruncCylinder'

      case(Z_AXIS)
        tst % axis = Z_AXIS
        tst % plane = [X_AXIS, Y_AXIS]
        type = 'zTruncCylinder'

      case default
        print *, "Should not happen. Wrong direction in testcase constructor"
        error stop

    end select

    ! Set origin & halfwidth
    origin = TWO
    origin(tst % plane(1)) = ONE
    origin(tst % plane(2)) = TWO

    hw = 1.5_defReal
    r = 2.0_defReal

    ! Build surface
    call dict % init(5)
    call dict % store('id', 75)
    call dict % store('type', type)
    call dict % store('origin', origin)
    call dict % store('halfwidth', hw)
    call dict % store('radius', r)
    call tst % surf % init(dict)

  end function newTestCase

  !!
  !! Deconstruct the test case
  !!
  subroutine tearDown(this)
    class(test_truncCylinder), intent(inout) :: this

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
    class(test_truncCylinder), intent(inout) :: this
    real(defReal), dimension(6) :: aabb, ref
    character(nameLen)          :: name
    integer(shortInt)                :: ax, p1, p2
    real(defReal), parameter    :: TOL = 1.0E-6_defReal

    ! Get axis and diffrent planar directions
    ax = this % axis
    p1 = this % plane(1)
    p2 = this % plane(2)

    ! Test ID
    @assertEqual(75, this % surf % id())

    ! Change ID
    call this % surf % setID(1)
    @assertEqual(1, this % surf % id())

    ! Name
    select case(this % axis)
      case(X_AXIS)
        name = 'xTruncCylinder'
      case(Y_AXIS)
        name = 'yTruncCylinder'
      case(Z_AXIS)
        name = 'zTruncCylinder'
    end select
    @assertEqual(name, this % surf % myType())

    ! Bounding Box
    ref([ax, p1, p2]) = [0.5_defReal, -1.0_defReal, 0.0_defReal]
    ref([ax, p1, p2]+3) = [3.5_defReal, 3.0_defReal, 4.0_defReal]
    aabb = this % surf % boundingBox()
    @assertEqual(ref, aabb, TOL)

  end subroutine testMisc

  !!
  !! Test boundary conditions
  !!
@Test(cases=[1,2,3])
  subroutine testBC(this)
    class(test_truncCylinder), intent(inout) :: this
    integer(shortInt)                :: ax, p1, p2
    integer(shortInt), dimension(3)  :: pe
    real(defReal), dimension(3)      :: r, u, r_ref, u_ref
    real(defReal), parameter :: TOL = 1.0E-6

    ! Get axis and diffrent planar directions
    ax = this % axis
    p1 = this % plane(1)
    p2 = this % plane(2)
    pe = [ax, p1, p2] ! Common permutation

    ! *** CASE 1 VACUUM + REFLECTIVE
    call this % surf % setBC([VACUUM_BC, REFLECTIVE_BC])

    ! ** Explicit BC
    ! Vacuum face
    r(pe) = [0.5_defReal, 0.0_defReal, 2.3_defReal]
    u(pe) = [-ONE, ZERO, ZERO]
    r_ref = r
    u_ref = u
    call this % surf % explicitBC(r, u)
    @assertEqual(r_ref, r, TOL)
    @assertEqual(u_ref, u, TOL)

    ! Reflective face
    r(pe) = [3.5_defReal, 0.0_defReal, 2.3_defReal]
    u(pe) = [ONE, ZERO, ZERO]
    r_ref = r
    u_ref = u
    u_ref(ax) = -u(ax)
    call this % surf % explicitBC(r, u)
    @assertEqual(r_ref, r, TOL)
    @assertEqual(u_ref, u, TOL)

    ! ** Transform BC
    r(pe) = [30.0_defReal, 0.0_defReal, 2.3_defReal]
    u(pe) = [ONE, ZERO, ZERO]
    r_ref(pe) = [-23.0_defReal, 0.0_defReal, 2.3_defReal]
    u_ref(pe) = [-ONE, ZERO, ZERO]
    call this % surf % transformBC(r, u)
    @assertEqual(r_ref, r, TOL)
    @assertEqual(u_ref, u, TOL)

    ! **** CASE 2 PERIODIC
    call this % surf % setBC([PERIODIC_BC, PERIODIC_BC])

    ! Explicit BC
    ! Periodic face
    r(pe) = [3.5_defReal, 0.0_defReal, 2.3_defReal]
    u(pe) = [ONE, ZERO, ZERO]
    r_ref(pe) = [0.5_defReal, 0.0_defReal, 2.3_defReal]
    u_ref = u
    call this % surf % explicitBC(r, u)
    @assertEqual(r_ref, r, TOL)
    @assertEqual(u_ref, u, TOL)

    ! ** Transform BC
    r(pe) = [12.0_defReal, 0.0_defReal, 2.3_defReal]
    u(pe) = [ONE, ZERO, ZERO]
    r_ref(pe) = [3.0_defReal, 0.0_defReal, 2.3_defReal]
    u_ref = u
    call this % surf % transformBC(r, u)
    @assertEqual(r_ref, r, TOL)
    @assertEqual(u_ref, u, TOL)

  end subroutine testBC

  !!
  !! Test Halfspaces membership
  !!
@Test(cases=[1,2,3])
  subroutine testHalfspace(this)
    class(test_truncCylinder), intent(inout) :: this
    integer(shortInt)                         :: ax, p1, p2
    integer(shortInt), dimension(3)           :: pe
    real(defReal), dimension(3)               :: r, u
    real(defReal)                             :: eps

    ! Get axis and diffrent planar directions
    ax = this % axis
    p1 = this % plane(1)
    p2 = this % plane(2)
    pe = [ax, p1, p2] ! Common permutation

    ! ** Outside
    r(pe) = [ 5.0_defReal, 3.0_defReal, 4.0_defReal ]
    u(pe) = [ZERO, ONE, ZERO]
    @assertTrue(this % surf % halfspace(r, u))

    ! Diffrent position
    r(pe) = [1.3_defReal, 1.0_defReal, 5.0_defReal]
    @assertTrue(this % surf % halfspace(r, u))

    ! One more
    r(pe) = [1.3_defReal, -2.0_defReal, 3.3_defReal]
    @assertTrue(this % surf % halfspace(r, u))

    ! ** Inside the surface
    ! Close to plane
    r(pe) = [3.3_defReal, 1.1_defReal, 2.2_defReal]
    @assertFalse(this % surf % halfspace(r, u))

    ! Close to cylinder
    r(pe) = [2.1_defReal, -0.7_defReal, 2.1_defReal]
    @assertFalse(this % surf % halfspace(r, u))

    ! ** Close to surface
    ! Cylinder
    r(pe) = [3.3_defReal, -1.0_defReal, 2.0_defReal]
    u(pe) = [ONE, ONE, ZERO]
    u = u/norm2(u)
    @assertFalse(this % surf % halfspace(r, u))
    @assertTrue(this % surf % halfspace(r, -u))

    ! Outside within surface tolerance
    eps = HALF * SURF_TOL
    r(pe) = [3.3_defReal, -1.0_defReal-eps, 2.0_defReal]
    u(pe) = [ONE, ONE, ZERO]
    u = u/norm2(u)
    @assertFalse(this % surf % halfspace(r, u))

    ! Outside outside surface tol
    eps = 1.001_defReal * SURF_TOL
    r(pe) = [3.3_defReal, -1.0_defReal-eps, 2.0_defReal]
    u(pe) = [ONE, ONE, ZERO]
    u = u/norm2(u)
    @assertTrue(this % surf % halfspace(r, u))

    ! Parallel within tolerance
    eps = HALF * SURF_TOL
    u(pe) = [ONE, ZERO, ZERO]

    r(pe) = [3.3_defReal, -1.0_defReal-eps, 2.0_defReal]
    @assertTrue(this % surf % halfspace(r, u))

    r(pe) = [3.3_defReal, -1.0_defReal+eps, 2.0_defReal]
    @assertFalse(this % surf % halfspace(r, u))

    ! Plane
    r(pe) = [3.5_defReal, 1.3_defReal, 1.8_defReal]
    u(pe) = [-ONE, ZERO, ONE]
    u = u/norm2(u)
    @assertFalse(this % surf % halfspace(r, u))
    @assertTrue(this % surf % halfspace(r, -u))

    ! Outside within surface tolerance
    eps = HALF * SURF_TOL
    r(pe) = [3.5_defReal+eps, 1.3_defReal, 1.8_defReal]
    u(pe) = [-ONE, ONE, ZERO]
    u = u/norm2(u)
    @assertFalse(this % surf % halfspace(r, u))

    ! Outside outside surface tol
    eps = 1.001_defReal * SURF_TOL
    r(pe) = [3.5_defReal+eps, 1.3_defReal, 1.8_defReal]
    u(pe) = [-ONE, ONE, ZERO]
    u = u/norm2(u)
    @assertTrue(this % surf % halfspace(r, u))

    ! Parallel within tolerance
    eps = HALF * SURF_TOL
    u(pe) = [ZERO, -ONE, ONE]
    u = u/norm2(u)

    r(pe) = [3.5_defReal+eps, 1.3_defReal, 1.8_defReal]
    @assertTrue(this % surf % halfspace(r, u))

    r(pe) = [3.5_defReal-eps, 1.3_defReal, 1.8_defReal]
    @assertFalse(this % surf % halfspace(r, u))

  end subroutine testHalfspace

  !!
  !! Test distance calculation
  !!
@Test(cases=[1, 2, 3])
  subroutine testDistance(this)
    class(test_truncCylinder), intent(inout) :: this
    integer(shortInt)                   :: ax, p1, p2
    integer(shortInt), dimension(3)     :: pe
    real(defReal), dimension(3)         :: r, u
    real(defReal)                       :: ref, eps
    real(defReal), parameter :: TOL = 1.0E-7

    ! Get axis and diffrent planar directions
    ax = this % axis
    p1 = this % plane(1)
    p2 = this % plane(2)
    pe = [ax, p1, p2] ! Common permutation

    ! ** Outside
    ! *Plane hits
    r(pe) = [5.0_defReal, 1.0_defReal, 2.0_defReal ]

    ! Direct hit
    u(pe) = [-ONE, ZERO, ZERO]
    ref = 1.5_defReal
    @assertEqual(ref, this % surf % distance(r, u), ref * TOL)

    ! Oblique hit
    u(pe) = [-ONE, -ONE, ZERO]
    u = u/norm2(u)
    ref = 1.5_defReal * SQRT2
    @assertEqual(ref, this % surf % distance(r, u), ref * TOL)

    ! Miss
    u(pe) = [-ONE, -TWO, -TWO]
    u = u/norm2(u)
    @assertEqual(INF, this % surf % distance(r, u))

    ! Parallel miss
    u(pe) = [ZERO, -TWO, -TWO]
    u = u/norm2(u)
    @assertEqual(INF, this % surf % distance(r, u))

    ! Wrong direction
    u(pe) = [ONE, -ONE, -ONE]
    u = u/norm2(u)
    @assertEqual(INF, this % surf % distance(r, u))

    ! * Cylinder hits
    r(pe) = [1.5_defReal, -2.0_defReal, 2.0_defReal ]

    ! Direct hit
    u(pe) = [ZERO, ONE, ZERO]
    ref = 1.0_defReal
    @assertEqual(ref, this % surf % distance(r, u), ref * TOL)

    ! Oblique hit 1
    u(pe) = [ONE, ONE, ZERO]
    u = u/norm2(u)
    ref = 1.0_defReal * SQRT2
    @assertEqual(ref, this % surf % distance(r, u), ref * TOL)

    ! Oblique hit 2
    u(pe) = [ZERO, 1.5_defReal, sqrt(7.0_defReal)* HALF]
    ref = norm2(u)
    u = u/norm2(u)
    @assertEqual(ref, this % surf % distance(r, u), ref * TOL)

    ! Miss
    u(pe) = [ONE, ONE, TWO]
    u = u/norm2(u)
    @assertEqual(INF, this % surf % distance(r, u))

    ! Parallel miss
    u(pe) = [ONE, ZERO, ZERO]
    @assertEqual(INF, this % surf % distance(r, u))

    ! Wrong direction
    u(pe) = [ONE, -ONE, ZERO]
    @assertEqual(INF, this % surf % distance(r, u))

    ! ** At the surface - Plane
    ! Outside within surf tolerance
    eps = HALF * SURF_TOL
    r(pe) = [3.5_defReal+eps, 1.3_defReal, 1.8_defReal]
    u(pe) = [-ONE, ZERO, ZERO]
    ref = 3.0_defReal
    @assertEqual(ref, this % surf % distance(r, u), ref * TOL)

    ! Inside going outside
    r(pe) = [3.5_defReal-eps, 1.3_defReal, 1.8_defReal]
    @assertEqual(INF, this % surf % distance(r, -u))

    ! ** At the surface - Cylinder
    ! Outside within surface tolerance
    eps = HALF * SURF_TOL
    r(pe) = [3.3_defReal, -1.0_defReal-eps, 2.0_defReal]
    u(pe) = [ZERO, ONE, ZERO]
    ref = 4.0_defReal
    @assertEqual(ref, this % surf % distance(r, u), ref * TOL)

    ! Inside going outside
    r(pe) = [3.3_defReal, -1.0_defReal+eps, 2.0_defReal]
    @assertEqual(INF, this % surf % distance(r, -u))

    ! ** Inside
    r(pe) = [3.0_defReal, 0.0_defReal, 2.0_defReal]

    ! Oblique hits
    u(pe) = [-ONE, -ONE, ZERO]
    u = u/norm2(u)
    ref = SQRT2
    @assertEqual(ref, this % surf % distance(r, u), ref * TOL)

    ref = HALF * SQRT2
    @assertEqual(ref, this % surf % distance(r, -u), ref * TOL)

    ! Parallel Planes
    u(pe) = [ZERO, -ONE, ZERO]
    ref = 1.0_defReal
    @assertEqual(ref, this % surf % distance(r, u), ref * TOL)

    ! Parallel Cylinder
    u(pe) = [-ONE, ZERO, ZERO]
    ref = 2.5_defReal
    @assertEqual(ref, this % surf % distance(r, u), ref * TOL)

  end subroutine testDistance

   !!
   !! Test Edge Cases
   !!
   !! Test unlikley cases to make sure that halfspace + distance
   !! procedures allow particle to escape
   !!
 @Test(cases=[1, 2, 3])
   subroutine testEdgeCases(this)
     class(test_truncCylinder), intent(inout) :: this
     real(defReal), dimension(3)    :: r, u
    integer(shortInt), dimension(3) :: pe
     integer(shortInt)              :: ax, p1, p2
     real(defReal)                  :: eps, d
     logical(defBool)               :: hs

     ! Get axis and diffrent planar directions
     ax = this % axis
     p1 = this % plane(1)
     p2 = this % plane(2)
     pe = [ax, p1, p2] ! Common permutation

     ! ** Try to escape the surface beeing in a corner
     ! * Particle is almost at the corner
     !   Either it is outside or can escape with a short movment in next step
     eps =  5.0_defReal * epsilon(eps)
     r([ax, p1, p2]) = [3.5_defReal-eps, -1.0_defReal+eps, 2.0_defReal]
     u([ax, p1, p2]) = [-TWO, -ONE, ZERO]
     u = u/norm2(u)
     hs = this % surf % halfspace(r, u)
     if (.not.hs) then ! Perform small movment
       d = this % surf % distance(r, u)
       @assertTrue( abs(d) < 1.0E-6)
       @assertTrue(this % surf % halfspace(r + d*u, u))
     end if

     ! Try other direction
     u([ax, p1, p2]) = [ONE, TWO, ZERO]
     u = u/norm2(u)
     hs = this % surf % halfspace(r, u)
     if (.not.hs) then ! Perform small movment
       d = this % surf % distance(r, u)
       @assertTrue( abs(d) < 1.0E-6)
       @assertTrue(this % surf % halfspace(r + d*u, u))
     end if

     ! * Asymetric corner
     eps =  5.0_defReal * epsilon(eps)
     r([ax, p1, p2]) = [3.5_defReal-eps, -1.0_defReal+TWO*eps, 2.0_defReal]
     u([ax, p1, p2]) = [-TWO, -ONE, ZERO]
     u = u/norm2(u)
     hs = this % surf % halfspace(r, u)
     if (.not.hs) then ! Perform small movment
       d = this % surf % distance(r, u)
       @assertTrue( abs(d) < 1.0E-6)
       @assertTrue(this % surf % halfspace(r + d*u, u))
     end if

     ! Try other direction
     r([ax, p1, p2]) = [3.5_defReal-TWO*eps, -1.0_defReal+eps, 2.0_defReal]
     u([ax, p1, p2]) = [ONE, TWO, ZERO]
     u = u/norm2(u)
     hs = this % surf % halfspace(r, u)
     if (.not.hs) then ! Perform small movment
       d = this % surf % distance(r, u)
       @assertTrue( abs(d) < 1.0E-6)
       @assertTrue(this % surf % halfspace(r + d*u, u))
     end if

   end subroutine testEdgeCases

end module truncCylinder_test
