program test
  
  use numPrecision
  use RNG_class, only : RNG


  implicit none
  type(RNG)         :: rand1
  integer(shortInt) :: i

  call rand1 % init(658758_longInt)

  print '(B64)', huge(1_8)
  print '(B64)', ibclr(huge(0_8),63)
!  do i=1,1000
!    print *, rand1 % get()
!
!  end do


end program test




