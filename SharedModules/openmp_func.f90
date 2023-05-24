!!
!! Contains wrappers around key OpenMP procedures
!!
!! Wrappers are required to enable calling "OpenMP" procedures in a code that may be compiled
!! without OpenMP support. The necessary conditional compilation is hidden from the user here
!!
!! Convention for wrappers is that the snake_case is changedto lowerCamelCase.
!!
!!
!!
module openmp_func
  use numPrecision
  use omp_lib
  implicit none
  private

  public :: ompGetThreadNum
  public :: ompGetMaxThreads
  public :: ompGetNumThreads
  public :: ompSetNumThreads

contains

  !!
  !! Get a number of a thread
  !!
  !! Result:
  !!   Number of a thread that is running [counting from 0!]
  !!
  function ompGetThreadNum() result(num)
    integer(shortInt) :: num
#ifdef _OPENMP
    num = omp_get_thread_num()
#else
    num = 0
#endif
  end function ompGetThreadNum

  !!
  !! Get the maximum number of threads of the program
  !!
  !! Result:
  !!   Maximum number of threads
  !!
  function ompGetMaxThreads() result(num)
    integer(shortInt) :: num
#ifdef _OPENMP
    num = omp_get_max_threads()
#else
    num = 1
#endif
  end function ompGetMaxThreads

  !!
  !! Return the number of active threads
  !!
  !! Result:
  !!  The number of threads currently active
  !!
  function ompGetNumThreads() result(num)
    integer(shortInt) :: num
#ifdef _OPENMP
    num = omp_get_num_threads()
#else
    num = 1
#endif
  end function ompGetNumThreads

  !!
  !! Set the maximum number of threads
  !!
  !! Does nothing if OpenMP is no enabled
  !!
  !! Args:
  !!   num [in] -> Maximum number of threads >= 1
  !!
  !! Error:
  !!   If num is <= 0. Sets maximum threads to 1
  !!
  subroutine ompSetNumThreads(num)
    integer(shortInt), intent(in) :: num
    integer(shortInt)             :: num_loc
#ifdef _OPENMP
    num_loc = max(1, num)
    call omp_set_num_threads(num_loc)
#else
    num_loc = 2 !Avoid compiler warning
#endif
  end subroutine ompSetNumThreads

end module openmp_func
