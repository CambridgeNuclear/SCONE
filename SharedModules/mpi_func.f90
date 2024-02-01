module mpi_func
  use numPrecision
#ifdef MPI
  use mpi_f08
#endif
  use genericProcedures, only : numToChar
  use errors_mod,        only : fatalError
  implicit none

  integer(shortInt), private :: worldSize
  integer(shortInt), private :: rank
  integer(shortInt), parameter  :: MASTER_RANK = 0

  !! Common MPI types
  type(MPI_Datatype)            :: MPI_DEFREAL

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

    call mpi_type_create_f90_real(precision(1.0_defReal), range(1.0_defReal), &
                                  MPI_DEFREAL, ierr)

    call mpi_type_commit(MPI_DEFREAL, ierr)

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

    call mpi_finalize(ierr)

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
  function isMPIMaster()
    logical(defBool) :: isMPIMaster

    isMPIMaster = (rank == 0)

  end function isMPIMaster

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
