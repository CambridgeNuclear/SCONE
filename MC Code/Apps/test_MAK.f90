program test
  
  use numPrecision
  use linearAlgebra_func


  implicit none
  real(8), dimension(2,2) :: A
  real(8), dimension(2),target :: s,f, x
!  real(8), dimension(:), pointer :: ptr
!  real(8), dimension(2)   :: WR, WI
!  real(8), dimension(2,2) :: VL, VR
!  real(8), dimension(10)  :: WORK
!  integer(shortInt) :: info
!

  ! Set linear system
  A(1,:) = [0.9_defReal,  0.1_defReal]
  A(2,:) = [0.1_defReal, 0.9_defReal]

  s = [ONE, -ONE] / sqrt(TWO)
  f = [ONE, ONE] / sqrt(TWO)

  x = ZERO
  call dgemv('N',2,2,1.0_defReal, A, 2, [ONE, ONE] , 1, 1.0_defReal, x, 1)
  print *, x

  call solveAdjointProblem(A,x,s,f)
  print *, x


!  print *, info
!  print *, A
!  print *, WR
!  print *, WI
!  !print *, VL
!  print *, VR



end program test




