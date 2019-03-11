!!
!! This module is a frontend for a linear algebra libraries
!! No libraries should be used by other modules directly. If an additional procedure is needed
!! approperiate interface should be added here.
!!
module linearAlgebra_func

  use numPrecision
  use genericProcedures, only : fatalError, numToChar
  use iso_fortran_env, only : real64, int32

  implicit none
  private

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Public Module interface
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  public :: eig
  public :: kill_linearAlgebra


!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! External LAPACK Procedures interfaces
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! lapack_geev
  !! LAPACK General Matrix Eigenvalue Solver (GE- General Matrix; EV - Eigenvalue)
  !! For full documentation of these procedures refer to:
  !! http://www.netlib.org/lapack/explore-html/index.html
  !! (Just search the webpage for dgeev or sgeev and you will find a very clear doc)
  !!
  !! Interface:
  !! jobvl [in]   -> 'V' or 'N'. With 'V' computes left eigenvectors
  !! jobvr [in]   -> 'V' or 'N'. With 'V' computes right eigenvectors
  !! N     [in]   -> Order of the matrix (NxN). N >=0
  !! A     [inout]-> LDA x N Matrix. Will be changed in the algorithm!
  !! LDA   [in]   -> Leading size of A. LDA >= max(1,N). Must be N
  !! WR    [out]  -> Vector of size N. Real Part of eigenvalues
  !! WI    [out]  -> Vector of size N. Imaginary parts of eigenvalues
  !! VL    [out]  -> LDVL x N Matrix of left eigenvalues. If JOBVL='N' it is not referenced
  !! LDVL  [in]   -> Leading dimension of VL
  !! VR    [out]  -> LDVR x N Matrix of left eigenvalues. If JOBVR='N' it is not referenced
  !! LDVR  [in]   -> Leading dimension of VR
  !! WORK  [out]  -> Work space. If LWORK = -1. On exit WORK(1) is size of optimal workspace
  !! LWORK [in]   -> Size of workspace. LWORK >= 3N for pure eigenvalue calculation LWORK >= 4N
  !!                 if any eigenvectors are also requested. If LWORK = -1 optimal size query.
  !! INFO  [out]  -> Error flag. INFO = 0 for succesfull exectution
  !!
  !! Note that the only difference in the dgeev and sgeev is kinf of the real arguments. Rest of
  !! the definition is identical.
  !!
  interface lapack_geev

    !!
    !! Double precision. 64-bit float
    !!
    subroutine dgeev(jobvl, jobvr, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO)
      use iso_fortran_env, only : real64, int32
      implicit none
      character(1), intent(in)                         :: jobvl
      character(1), intent(in)                         :: jobvr
      integer(int32),intent(in)                        :: N
      integer(int32), intent(in)                       :: LDA
      real(real64),dimension(LDA,N),intent(inout)      :: A
      real(real64),dimension(N), intent(out)           :: WR
      real(real64),dimension(N), intent(out)           :: WI
      integer(int32), intent(in)                       :: LDVL
      real(real64),dimension(LDVL,N), intent(out)      :: VL
      integer(int32),intent(in)                        :: LDVR
      real(real64),dimension(LDVR,N),intent(out)       :: VR
      integer(int32),intent(in)                        :: LWORK
      real(real64),dimension(max(1,LWORK)),intent(out) :: WORK
      integer(int32),intent(out)                       :: INFO
    end subroutine dgeev

    !!
    !! Single precision. 32-bit float
    !!
    subroutine sgeev(jobvl, jobvr, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO)
      use iso_fortran_env, only : real32, int32
      implicit none
      character(1), intent(in)                         :: jobvl
      character(1), intent(in)                         :: jobvr
      integer(int32),intent(in)                        :: N
      integer(int32), intent(in)                       :: LDA
      real(real32),dimension(LDA,N),intent(inout)      :: A
      real(real32),dimension(N), intent(out)           :: WR
      real(real32),dimension(N), intent(out)           :: WI
      integer(int32), intent(in)                       :: LDVL
      real(real32),dimension(LDVL,N), intent(out)      :: VL
      integer(int32),intent(in)                        :: LDVR
      real(real32),dimension(LDVR,N),intent(out)       :: VR
      integer(int32),intent(in)                        :: LWORK
      real(real32),dimension(max(1,LWORK)),intent(out) :: WORK
      integer(int32),intent(out)                       :: INFO
    end subroutine sgeev
  end interface lapack_geev

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Module VARIABLES
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !! This variable must be private to each OpenMP thread
  real(defReal),dimension(:),allocatable,target :: workspace

contains

  !!
  !! Calculates real part of eigenvalues and
  !! right eigenvectors of a square matrix A
  !!
  !! A - any square real matrix
  !! k - vector of real part of the eigenvalue
  !! V - matrix of right eigenvectors
  !!
  !! For +ve determinant matrix eigenvalues are in descending order
  !! For -ve determinant matrix eigenvalues are in ascending order
  !! Eigenvectors are normalised to Euclidean norm of 1.0
  !!
  subroutine eig(k, V, A)
    real(defReal),dimension(:),intent(out)    :: k
    real(defReal),dimension(:,:), intent(out) :: V
    real(defReal),dimension(:,:), intent(in)  :: A
    integer(int32)                            :: N, mem, st, info
    real(defReal),dimension(:,:),pointer      :: A_t, VR_t, VL_t
    real(defReal),dimension(:),pointer        :: Re, Im, Work
    character(100),parameter :: Here = 'eig (linearAlgebra_func.f90)'

    ! Verify size of inputs
    N = size(A,1)
    if( any(shape(V) /= N)) then
      call fatalError(Here,'Invalid shape of eigenvector result array. Is not NxN')

    else if(size(k) /= N) then
      call fatalError(Here,'Invalid size of eigenvalue result vectorr. Is not size N')

    else if ( any(shape(A) /= N)) then
      call fatalError(Here,'Invalid shape of array A. Is not NxN')

    end if

    ! Calculate memory required and ensure that memory is avalible
    ! Mem for: A     V    VL   k    Work
    mem    =  N*N + N*N + N + 2*N + 5*N

    ! Ensure that memory is avalible
    call getMem(mem)

    ! Associate workspace memory with different variables
    ! Use pointers to change ranks
    st = 1
    A_t(1:N,1:N)  => workspace(st : st + N*N-1)

    st = st + N*N
    VR_t(1:N,1:N) => workspace(st : st + N*N-1)

    st = st + N*N
    VL_t(1:1,1:N) => workspace(st : st + N-1)

    st = st + N
    Re   => workspace(st : st + N-1)

    st = st + N
    Im   => workspace(st : st + N-1)

    st = st + N
    Work => workspace(st : size(workspace))

    ! Copy Input on the workmemory
    A_t = A

    ! Perform calculation
    call lapack_geev('N','V', N, A_t, N, Re, Im, VL_t, 1, VR_t, N, Work, size(Work), info)

    if( info /= 0) then
      call fatalError(Here,'LINPACK procedure failed with error: '//numToChar(info))
    end if

    ! Copy the results out
    V = VR_t
    k = Re

  end subroutine eig


  !!
  !! Makes sure that workspace is allocated and has size > N
  !!
  subroutine getMem(N)
    integer(shortInt), intent(in) :: N

    if(allocated(workspace)) then
      ! Check that workspace has sufficient size
      if(size(workspace) < N) then
        deallocate(workspace)
        allocate(workspace(N))
      end if

    else
      allocate(workspace(N))

    end if

  end subroutine getMem


  !!
  !! Returns module to its uninitialised state
  !!
  subroutine kill_linearAlgebra()

    if(allocated(workspace)) deallocate(workspace)

  end subroutine kill_linearAlgebra

end module linearAlgebra_func
