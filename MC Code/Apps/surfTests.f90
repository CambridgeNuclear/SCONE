program surfTests
  implicit none

  call xSquareCylTests()
  !call boxTest()

end program surfTests


subroutine boxTest()
  use numPrecision
  use universalVariables
  use genericProcedures, only : rotateVector
  use vector_class,      only : vector
  use box_class,         only : box

  implicit none

  type(vector) :: r
  type(vector) :: u
  type(vector) :: vt
  type(box)    :: tBox
  real(defReal) :: ev, mu, phi
  real(defReal) :: temp
  integer(shortInt) :: i

  call tBox % init([ZERO, ZERO, ZERO], [0.5_8, 0.6_8, 0.7_8],1)
  !call tBox % setBoundaryConditions([1,1,2,1,1,0])

  temp = -1.0_8 * surface_tol - 0.7_8
  r = [0.3_8 ,-0.6_8 ,temp]
  call tBox % evaluate(ev,r)
  print *, ev

  temp = -1.1_8 * surface_tol - 0.7_8
  u =[ ZERO, ZERO, ONE]
  r = [0.3_8 ,-0.6_8 ,temp]
  print *, "+VE HALFSPACE: ", tBox % halfspace(r,u)

  temp =  -0.7_8
  u =[ ZERO, ZERO, ONE]
  r = [0.3_8 ,-0.2_8 ,temp]
  print *, "-VE HALFSPACE: ", tBox % halfspace(r,u)

  temp = - 0.7_8
  u =[ ZERO, ZERO, -ONE]
  r = [0.3_8 ,-0.2_8 ,temp]
  print *, "+VE HALFSPACE: ", tBox % halfspace(r,u)
  vt = tBox % normalVector(r)
  print *, vt % v


  print *, "DISTANCE CHECKS"
  mu = ONE-0.00000001_8
  !mu = 0.0_8

  u = [mu, -(ONE-mu*mu), 0.0000000_8]
  u = u / u % L2norm()

  r = [ 0.0_8 -0.0_8 * surface_tol ,&
        0.6_8 +1.005_8 * surface_tol ,&
        0.0_8 -0.0_8 * surface_tol ]

  print *, r
  print *, "HALFSPACE: ", tBox % halfspace(r,u)
  call tBox % distance(ev,i,r,u)
  print *, "DISTANCE : ", ev, i


  print *, "BC CHECK"

  call tBox % init([0.5_8, ZERO, ZERO], [0.5_8, 0.6_8, 0.7_8],1)
  call tBox % setBoundaryConditions([1,1,1,1,2,2])

  r = [ -0.5_8   -0.0_8 * surface_tol ,&
        -0.6_8   +0.0_8 * surface_tol ,&
         0.7_8   -0.0_8 * surface_tol ]

  u = [mu, -(ONE-mu*mu), 1.0000000_8]
  u = u / u % L2norm()

  print *, "BEFORE"
  print *, "r: ", r % v
  print *, "u: ", u % v


  call tBox % boundaryTransform(r,u)
  print *, "AFTER"
  print *, "r: ", r % v
  print *, "u: ", u % v

end subroutine  boxTest


subroutine xSquareCylTests()
  use numPrecision
  use universalVariables
  use genericProcedures,        only : rotateVector
  use vector_class,             only : vector
  use xSquareCylinder_class,    only : xSquareCylinder

  implicit none

  type(vector)           :: r
  type(vector)           :: u
  type(vector)           :: vt
  type(xSquareCylinder)  :: tSq
  real(defReal)          :: ev, mu, phi
  real(defReal)          :: temp
  integer(shortInt)      :: i
  character(:), allocatable :: defString
  integer(shortInt), dimension(6) :: BC


  call tSq % init( [0.0_8, -0.2_8, 0.1_8], [0.0_8, 0.6_8, 0.5_8], 7)
  call tSq % getDef(defString)
  print *, repeat('<>',50)
  print *, defString

  print *, "***HALFSPACE CHECKS***"
  print *, "CASE INSIDE"
  mu = -0.5
  r = [17.1_8, 0.2_8, 0.1_8]
  u = [mu,(1-mu*mu), 0.3_8]
  u = u/u % L2norm()
  print *, "r: ", r % v
  print *, "u: ", u % v
  print *, "-VE HALFSPACE: ", tSq % halfspace(r,u)

  print *, "CASE AT THE SURFACE -> GOING INSIDE"
  mu = -0.5
  r = [17.1_8, -0.8_8, 0.1_8]
  u = [mu,(1-mu*mu), 0.3_8]
  u = u/u % L2norm()
  print *, "r: ", r % v
  print *, "u: ", u % v
  print *, "-VE HALFSPACE: ", tSq % halfspace(r,u)

  print *, "CASE AT THE SURFACE -> GOING OUTSIDE"
  mu = -0.5
  r = [17.1_8, -0.8_8, 0.1_8]
  u = [mu,-(1-mu*mu), 0.3_8]
  u = u/u % L2norm()
  print *, "r: ", r % v
  print *, "u: ", u % v
  print *, "+VE HALFSPACE: ", tSq % halfspace(r,u)

  print *, "CASE AT THE SURFACE CLOSE TO CORNER SHALLOW ANGLE -> GOING INSIDE"
  mu = -0.5
  r = [17.1_8, -0.8_8, 0.6_8-1.1*surface_tol]
  u = [mu, 700.0_8, -0.3_8]
  u = u/u % L2norm()
  print *, "r: ", r % v
  print *, "u: ", u % v
  print *, "-VE HALFSPACE: ", tSq % halfspace(r,u)
  vt = tSq % normalVector(r)
  print *, vt %v

  print *, "CASE AT THE CORNER -> GOING INSIDE"
  mu = tiny(mu)
  r = [17.1_8, -0.8_8, 0.6_8]
  u = [mu, (1-mu*mu), -0.3_8]
  u = u/u % L2norm()
  print *, "r: ", r % v
  print *, "u: ", u % v
  print *, "-VE HALFSPACE: ", tSq % halfspace(r,u)
  vt = tSq % normalVector(r)
  print *, vt %v

  print *, "CASE OUTSIDE"
  r = [17.1_8, -0.8_8 +1.1_8*surface_tol, 0.6+1.1_8*surface_tol]
  u = [mu, (1-mu*mu), -0.3_8]
  u = u/u % L2norm()
  print *, "r: ", r % v
  print *, "u: ", u % v
  print *, "+VE HALFSPACE: ", tSq % halfspace(r,u)

  print *, "***DISTANCE CHECKS***"
  print *, "CASE INSIDE GOING INTO CORNER"
  mu = sqrt(TWO)*0.5
  r = [ -0.1_8, -0.1_8, 0.1_8]
  u = [0.0_8, mu, sqrt(1-mu*mu)]
  u = u/u % L2norm()

  print *, "r: ", r % v
  print *, "u: ", u % v
  call tSq % distance(ev,i,r,u)
  print *, "DISTANCE : ", ev, "SHOULD BE :", sqrt(TWO)*0.5, " AT SURFACE 3/5", i

  print *, "CASE AT THE SURFACE GOING VERY SHALLOW OUTSIDE "
  mu = 0.999_8
  r = [ -0.1_8, 0.0_8, 0.6_8+0.8_8*surface_tol]
  u = [0.0_8, mu, sqrt(1-mu*mu)]
  u = u/u % L2norm()

  print *, "r: ", r % v
  print *, "u: ", u % v
  call tSq % distance(ev,i,r,u)
  print *, "DISTANCE : ", ev, "SHOULD BE :", INFINITY, " AT SURFACE 0", i

  print *, "CASE OUTSIDE MISS"
  mu = 0.999_8
  r = [ -0.1_8, 0.5_8, 0.7_8]
  u = [0.0_8, mu, sqrt(1-mu*mu)]
  u = u/u % L2norm()

  print *, "r: ", r % v
  print *, "u: ", u % v
  call tSq % distance(ev,i,r,u)
  print *, "DISTANCE : ", ev, "SHOULD BE :", INFINITY, " AT SURFACE 0", i

  print *, "CASE OUTSIDE HIT"
  mu = 0.999_8
  r = [ -0.1_8, -1.0_8, 0.7_8]
  u = [ZERO, ONE, ZERO]
  u = u/u % L2norm()

  print *, "r: ", r % v
  print *, "u: ", u % v
  call tSq % distance(ev,i,r,u)
  print *, "DISTANCE : ", ev, "SHOULD BE :", .2_8, " AT SURFACE 4", 4

  print *, "***BC CHECKS***"
  BC = [ 0, 0, 1, 0, 2, 2]
  call tSq % setBoundaryConditions(BC)
  print *, "BCs: ", BC

  print *, "REFLECTION FROM THE SURFACE"
  mu = 0.999_8
  r = [ -0.1_8, 0.4_8 + 0.9_8*surface_tol, 0.0_8]
  u = [ONE, ONE, ONE]
  u = u/u % L2norm()

  print *, "BEFORE"
  print *, "r: ", r % v
  print *, "u: ", u % v

  call tSq % boundaryTransform(r,u)
  print *, "AFTER"
  print *, "r: ", r % v
  print *, "u: ", u % v

  print *, " LEAK FROM THE SURFACE"
  mu = 0.999_8
  r = [ -0.1_8, -4.4_8, 0.0_8]
  u = [ONE, ONE, ONE]
  u = u/u % L2norm()

  print *, "BEFORE"
  print *, "r: ", r % v
  print *, "u: ", u % v

  call tSq % boundaryTransform(r,u)
  print *, "AFTER"
  print *, "r: ", r % v
  print *, "u: ", u % v

  print *, " TRANSLATION ACROSS SURFACES WITH REFLECTION"
  mu = 0.999_8
  r = [ -0.1_8, 0.5_8, 0.7_8]
  u = [ONE, ONE, ONE]
  u = u/u % L2norm()

  print *, "BEFORE"
  print *, "r: ", r % v
  print *, "u: ", u % v

  call tSq % boundaryTransform(r,u)
  print *, "AFTER"
  print *, "r: ", r % v
  print *, "u: ", u % v


  print *, repeat('<>',50)

end subroutine
