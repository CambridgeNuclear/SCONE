module mpi_func
  use numPrecision
#ifdef MPI
  use mpi_f08
#endif
  implicit none

  integer(shortInt), private :: worldSize
  integer(shortInt), private :: rank

contains

  !!
  !! Initialise MPI environment
  !!
  !! Needs to be called at the beginning of calculation before any MPI calls
  !!
  subroutine mpiInit()
#ifdef MPI
    integer(shortInt) :: ierr
    call mpi_init(ierr)

    call mpi_comm_size(MPI_COMM_WORLD, worldSize, ierr)

    call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)

#else
    worldSize = 1
    rank = 0
#endif
  end subroutine mpiInit

  !!
  !! Finalise MPI environment
  !!
  !! Needs to be called at the end of calculation after all MPI calls
  !!
  subroutine mpiFinalise()
#ifdef MPI
    integer(shortInt) :: ierr
    call MPI_Finalize(ierr)
#endif
  end subroutine mpiFinalise

  !!
  !! Get MPI world size
  !!
  !! It is the number of processes launched concurrently and communicating
  !! with each other
  !!
  function getMPIWorldSize() result(size)
    integer(shortInt) :: size
    size = worldSize
  end function getMPIWorldSize

  !!
  !! Return true if the process is the master process
  !!
  !! The master process is the one with rank 0
  !!
  function isMaster()
    logical(defBool) :: isMaster

    isMaster = (rank == 0)

  end function isMaster

  !!
  !! Get MPI rank
  !!
  !! It is the number of the process in the MPI world.
  !! Unlike in Fortran, it starts from 0
  !!
  function getMPIRank() result(r)
    integer(shortInt) :: r
    r = rank
  end function getMPIRank

end module mpi_func
