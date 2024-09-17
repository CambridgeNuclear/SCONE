module mpi_func
  use numPrecision
#ifdef MPI
  use mpi_f08
#endif
  use errors_mod,        only : fatalError

  implicit none

  integer(shortInt), private   :: worldSize = 1
  integer(shortInt), private   :: rank = 0
  integer(shortInt), parameter :: MASTER_RANK = 0

  !! Public type that replicates exactly particleState
  !!
  !! It is necessary for load balancing in the dungeon: particles have to be
  !! transferred betwen processes, and MPI doesn't allow to transfer types with
  !! type-bound procedures
  type, public :: particleStateDummy
    real(defReal)              :: wgt
    real(defReal),dimension(3) :: r
    real(defReal),dimension(3) :: dir
    real(defReal)              :: E
    integer(shortInt)          :: G
    logical(defBool)           :: isMG
    integer(shortInt)          :: type
    real(defReal)              :: time
    integer(shortInt)          :: matIdx
    integer(shortInt)          :: cellIdx
    integer(shortInt)          :: uniqueID
    integer(shortInt)          :: collisionN
    integer(shortInt)          :: broodID
  end type particleStateDummy

  !! Common MPI types
#ifdef MPI
  type(MPI_Datatype)   :: MPI_DEFREAL
  type(MPI_Datatype)   :: MPI_SHORTINT
  type(MPI_Datatype)   :: MPI_LONGINT
  type(MPI_Datatype)   :: MPI_PARTICLE_STATE
#endif

contains

  !!
  !! Initialise MPI environment
  !!
  !! Needs to be called at the beginning of calculation before any MPI calls
  !!
  subroutine mpiInit()
#ifdef MPI
    integer(shortInt)        :: ierr, stateSize
    type(particleStateDummy) :: state
    integer(kind = MPI_ADDRESS_KIND), dimension(:), allocatable :: displacements
    integer(shortInt), dimension(:), allocatable                :: blockLengths
    type(MPI_Datatype), dimension(:), allocatable               :: types

    call mpi_init(ierr)

    ! Read number of processes and rank of each process
    call mpi_comm_size(MPI_COMM_WORLD, worldSize, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)

    ! Define MPI type for DEFREAL
    call mpi_type_create_f90_real(precision(1.0_defReal), range(1.0_defReal), &
                                  MPI_DEFREAL, ierr)
    call mpi_type_commit(MPI_DEFREAL, ierr)

    ! Define MPI type for SHORTINT
    call mpi_type_create_f90_integer(range(1_shortInt), MPI_SHORTINT, ierr)
    call mpi_type_commit(MPI_SHORTINT, ierr)

    ! Define MPI type for LONGINT
    call mpi_type_create_f90_integer(range(1_longInt), MPI_LONGINT, ierr)
    call mpi_type_commit(MPI_LONGINT, ierr)

    ! Define MPI type for particleState
    ! Note that particleState has stateSize = 13 attributes; if an attribute is
    ! added to particleState, it had to be added here too
    stateSize = 13
    allocate(displacements(stateSize), blockLengths(stateSize), types(stateSize))

    ! Create arrays with dimension and type of each property of particleStateDummy
    blockLengths = (/1, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1/)
    types = (/MPI_DEFREAL, MPI_DEFREAL, MPI_DEFREAL, MPI_DEFREAL, MPI_SHORTINT, MPI_LOGICAL, MPI_SHORTINT, &
              MPI_DEFREAL, MPI_SHORTINT, MPI_SHORTINT, MPI_SHORTINT, MPI_SHORTINT, MPI_SHORTINT/)

    ! Create array of memory byte displacements
    call mpi_get_address(state % wgt, displacements(1), ierr)
    call mpi_get_address(state % r, displacements(2), ierr)
    call mpi_get_address(state % dir, displacements(3), ierr)
    call mpi_get_address(state % E, displacements(4), ierr)
    call mpi_get_address(state % G, displacements(5), ierr)
    call mpi_get_address(state % isMG, displacements(6), ierr)
    call mpi_get_address(state % type, displacements(7), ierr)
    call mpi_get_address(state % time, displacements(8), ierr)
    call mpi_get_address(state % matIdx, displacements(9), ierr)
    call mpi_get_address(state % cellIdx, displacements(10), ierr)
    call mpi_get_address(state % uniqueID, displacements(11), ierr)
    call mpi_get_address(state % collisionN, displacements(12), ierr)
    call mpi_get_address(state % broodID, displacements(13), ierr)
    displacements = displacements - displacements(1)

    ! Define new type
    call mpi_type_create_struct(stateSize, blockLengths, displacements, types, MPI_PARTICLE_STATE, ierr)
    call mpi_type_commit(MPI_PARTICLE_STATE, ierr)

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
  !! Get the share of work N for the current process
  !!
  !! Args:
  !!  N [in] -> Total measure of work (e.g. number of particles)
  !!
  !! Result:
  !!  The share of work for the current process
  !!
  function getWorkshare(N) result(share)
    integer(shortInt), intent(in) :: N
    integer(shortInt)             :: share

    share = (N + rank) / worldSize

  end function getWorkshare

  !!
  !! Get starting work offset for the current process
  !!
  !! Args:
  !!  N [in] -> Total measure of work (e.g. number of particles)
  !!
  !! Result:
  !!  The starting offset for the current process: offset = Sum_{i=0}^{rank-1} N_i
  !!  where N_i is the share of work for process i
  !!
  function getOffset(N) result(offset)
    integer(shortInt), intent(in) :: N
    integer(shortInt)             :: offset
    integer(shortInt)             :: remainder

    remainder = mod(N, worldSize)
    offset = N / worldSize * rank + max(0, remainder + rank - worldSize)

  end function getOffset

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
