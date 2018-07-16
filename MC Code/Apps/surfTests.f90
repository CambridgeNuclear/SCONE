program surfTests
  implicit none

  call boxTest()

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

  call tBox % init([ZERO, ZERO, ZERO], [0.5_8, 0.6_8, 0.7_8],1)
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
