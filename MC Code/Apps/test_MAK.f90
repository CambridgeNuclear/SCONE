program test
  
  use numPrecision
  use genericProcedures
  use timer_mod


  implicit none
  type(stopWatch) :: watch
  real(defReal), dimension(:),allocatable :: rand

  call watch % start()


  allocate(rand(100000000))


  call random_number(rand)

  call quickSort(rand)


  call watch % stop()

  print *, secToChar( watch % toSec())

end program test




