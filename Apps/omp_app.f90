program omp_app
  use numPrecision
  use openmp_func

  implicit none


  call ompSetNumThreads(2)

  !$omp parallel
  if(ompGetThreadNum() == 0) then
    print *, "Master thread. # of threads active ", ompGetNumThreads(), "of", ompGetMaxThreads()
  end if

  print *, 'Hello from thread: ', ompGetThreadNum()
  !$omp end parallel


end program omp_app
