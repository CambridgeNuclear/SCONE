!!
!! This module is a frontend for a linear algebra libraries
!! No libraries should be used by other modules directly. If an additional procedure is needed
!! approperiate interface should be added here.
!!
module linearAlgebra_func

  use numPrecision
  use genericProcedures, only : fatalError, numToChar
  use iso_fortran_env,   only : real64, int32

  implicit none
  private

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Public Module interface
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  public :: eig
  public :: solve
  public :: solveAdjointProblem
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
  !! Note that the only difference in the dgeev and sgeev is kind of the real arguments. Rest of
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

  !!
  !! lapack_gesv
  !! LAPACK General Matrix Linear Equation solver
  !! For full documentation of these procedures refer to:
  !! http://www.netlib.org/lapack/explore-html/index.html
  !! (Just search the webpage for dgesv or sgesv and you will find a very clear doc)
  !!
  !! Solves AX=B, where each column of X and VB corresponds to Ax=b.
  !! In outher words can solve for a number of vectors b simultaneously
  !! Uses LU decomposition with Partial Pivoting
  !!
  !! Interface:
  !! N    [in]    -> Number of linear equations to solve
  !! NRHS [in]    -> Number of "Right-hand Sides"
  !! A    [inout] -> On entry coefficient matrix A. On exit L and U factorisation: A=P*L*U
  !! LDA  [in]    -> Leading size of A. LDA >= max(1,N). Must be N
  !! IPIV [out]   -> Integer array of size(N). Pivot indices that define permulation matrix P;
  !!                 row i of the matrix was interchanged with row IPIV(i).
  !! B    [inout] -> Real LDB x NRHS array. On entry matrix B. On exit matrix X.
  !! LDB  [in]    -> Leading dimension of matrix B. Must be N.
  !! INFO [out]   -> Error flag. INFO = 0 for succesfull exectution
  !!
  !! Note that the only difference in the dgesv and sgesv is kind of the real arguments. Rest of
  !! the definition is identical.
  !!
  interface lapack_gesv
    !!
    !! Double precision. 64-bit float
    !!
    subroutine dgesv(N, NRHS, A, LDA, IPIV, B, LDB, INFO)
      use iso_fortran_env, only : real64, int32
      implicit none
      integer(int32), intent(in)                      :: N
      integer(int32), intent(in)                      :: NRHS
      integer(int32), intent(in)                      :: LDA
      real(real64),dimension(LDA,N), intent(inout)    :: A
      integer(int32),dimension(:), intent(out)        :: IPIV
      integer(int32), intent(in)                      :: LDB
      real(real64),dimension(LDB,NRHS), intent(inout) :: B
      integer(int32), intent(out)                     :: INFO
    end subroutine dgesv

    !!
    !! Single precision. 32-bit float
    !!
    subroutine sgesv(N, NRHS, A, LDA, IPIV, B, LDB, INFO)
      use iso_fortran_env, only : real32, int32
      implicit none
      integer(int32), intent(in)                      :: N
      integer(int32), intent(in)                      :: NRHS
      integer(int32), intent(in)                      :: LDA
      real(real32),dimension(LDA,N), intent(inout)    :: A
      integer(int32),dimension(:), intent(out)        :: IPIV
      integer(int32), intent(in)                      :: LDB
      real(real32),dimension(LDB,NRHS), intent(inout) :: B
      integer(int32), intent(out)                     :: INFO
    end subroutine sgesv
  end interface lapack_gesv

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! External BLAS Procedures interfaces
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! blas_dot
  !! BLAS Dot product procedure
  !! Note that BLAS interface is included in LAPACK documentation
  !! For full documentation of these procedures refer to:
  !! http://www.netlib.org/lapack/explore-html/index.html
  !! (Just search the webpage for dgesv or sgesv and you will find a very clear doc)
  !!
  !! Interface:
  !!   N [in]   -> Number of elements in vector inputs
  !!   DX [in]  -> First array of size (1 + (N-1) * abs(INCX))
  !!   INCX [in]-> Storage spacing between elements of DX
  !!   DY [in]  -> Second array of size (1 + (N-1) * abs(INCX))
  !!   INCY [in]-> Storage spacing between elements of DY
  !!
  !! Result:
  !!   Value of dot product in approperiate precision
  !!
  !! NOTE: Despite the fact it is not explicitly stated in the docs it is prudent to presume
  !!       that effective size of DX and DY, taking storage increpents under account, must be equal
  !!
  interface blas_dot
    !!
    !! Double precision. 64-bit float
    !!
    function ddot(N, DX, INCX, DY, INCY) result(DOT)
      use iso_fortran_env, only : real64, int32
      implicit none
      integer(int32), intent(in)           :: N
      real(real64),dimension(:),intent(in) :: DX
      integer(int32), intent(in)           :: INCX
      real(real64),dimension(:),intent(in) :: DY
      integer(int32), intent(in)           :: INCY
      real(real64)                         :: DOT
    end function ddot

    !!
    !! Single precision. 32-bit float
    !!
    function sdot(N, DX, INCX, DY, INCY) result(DOT)
      use iso_fortran_env, only : real32, int32
      implicit none
      integer(int32), intent(in)           :: N
      real(real32),dimension(:),intent(in) :: DX
      integer(int32), intent(in)           :: INCX
      real(real32),dimension(:),intent(in) :: DY
      integer(int32), intent(in)           :: INCY
      real(real32)                         :: DOT
    end function sdot
  end interface blas_dot

  !!
  !! blas_gemv
  !! BLAS General Matrix-Vector Multiplication
  !! Note that BLAS interface is included in LAPACK documentation
  !! For full documentation of these procedures refer to:
  !! http://www.netlib.org/lapack/explore-html/index.html
  !! (Just search the webpage for dgesv or sgesv and you will find a very clear doc)
  !!
  !! Computes one of the two:
  !!  1) y = alpha*A*x + beta*y
  !!  2) y = alpha*A**T*x + beta*y
  !!
  !! Interface:
  !!   TRANS [in]    -> character for 'N' computes 1). For 'T' computes 2)
  !!   M     [in]    -> Integer. Number of rows of matrix A. M >= 0
  !!   N     [in]    -> Integer. Number of columns of matrix A. N >= 0
  !!   ALPHA [in]    -> Real. Alpha in 1) and 2)
  !!   A     [in]    -> LDA x N Matrix
  !!   LDA   [in]    -> Leading dimension of matrix A. LDA >= max(1,M)
  !!   X     [in]    -> Vector x of size that matches N or M depending on TRANS
  !!   INCX  [in]    -> Increment for the elements of X
  !!   BETA  [in]    -> Real. Beta in 1) and 2)
  !!   Y     [inout] -> Vector y of size that matches N or M depending on TRANS
  !!   INCY  [in]    -> Increment for the elements of Y
  !!
  interface blas_gemv
    !!
    !! Double precision. 64-bit float
    !!
    subroutine dgemv(TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
      use iso_fortran_env, only : real64, int32
      implicit none
      character(1), intent(in)                    :: TRANS
      integer(int32), intent(in)                  :: M
      integer(int32), intent(in)                  :: N
      real(real64), intent(in)                    :: ALPHA
      integer(int32), intent(in)                  :: LDA
      real(real64), dimension(LDA, N), intent(in) :: A
      real(real64), dimension(:), intent(in)      :: X
      integer(int32), intent(in)                  :: INCX
      real(real64), intent(in)                    :: BETA
      real(real64), dimension(:), intent(inout)   :: Y
      integer(int32), intent(in)                  :: INCY
    end subroutine dgemv

    !!
    !! Singe precision. 32-bit float
    !!
    subroutine sgemv(TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
      use iso_fortran_env, only : real32, int32
      implicit none
      character(1), intent(in)                    :: TRANS
      integer(int32), intent(in)                  :: M
      integer(int32), intent(in)                  :: N
      real(real32), intent(in)                    :: ALPHA
      integer(int32), intent(in)                  :: LDA
      real(real32), dimension(LDA, N), intent(in) :: A
      real(real32), dimension(:), intent(in)      :: X
      integer(int32), intent(in)                  :: INCX
      real(real32), intent(in)                    :: BETA
      real(real32), dimension(:), intent(inout)   :: Y
      integer(int32), intent(in)                  :: INCY
    end subroutine sgemv
  end interface blas_gemv

  !!
  !! blas_axpy
  !! BLAS Constant * Vector + Vector Procedure
  !! Note that BLAS interface is included in LAPACK documentation
  !! For full documentation of these procedures refer to:
  !! http://www.netlib.org/lapack/explore-html/index.html
  !! (Just search the webpage for dgesv or sgesv and you will find a very clear doc)
  !!
  !! Computes:
  !!   y = DA * x + y
  !!
  !! Interface:
  !!   N    [in]    -> Size of vectors
  !!   DA   [in]    -> Real. Scalar DA.
  !!   DX   [in]    -> Vector x. Size > ( 1 + ( N - 1 )*abs( INCX )
  !!   INCX [in]    -> Storage spacing between elements of DX
  !!   DY   [inout] -> Vector y. Size > ( 1 + ( N - 1 )*abs( INCY )
  !!   INCY [in]    -> Storage spacing between elements of DY
  !!
  !! NOTE: It should be assumed that effective size of DX and DY, taking increments under account,
  !!       has to match.
  !!
  interface blas_axpy
    !!
    !! Double precision. 64-bit float
    !!
    subroutine daxpy(N, DA, DX, INCX, DY, INCY)
      use iso_fortran_env, only : real64, int32
      implicit none
      integer(int32), intent(in)                :: N
      real(real64), intent(in)                  :: DA
      real(real64), dimension(:), intent(in)    :: DX
      integer(int32), intent(in)                :: INCX
      real(real64), dimension(:), intent(inout) :: DY
      integer(int32), intent(in)                :: INCY
    end subroutine daxpy

    !!
    !! Single precision. 32-bit float
    !!
    subroutine saxpy(N, SA, SX, INCX, SY, INCY)
      use iso_fortran_env, only : real32, int32
      implicit none
      integer(int32), intent(in)                :: N
      real(real32), intent(in)                  :: SA
      real(real32), dimension(:), intent(in)    :: SX
      integer(int32), intent(in)                :: INCX
      real(real32), dimension(:), intent(inout) :: SY
      integer(int32), intent(in)                :: INCY
    end subroutine saxpy

  end interface blas_axpy

  !!
  !! blas_copy
  !! BLAS Procedure to copy vectors
  !! Note that BLAS interface is included in LAPACK documentation
  !! For full documentation of these procedures refer to:
  !! http://www.netlib.org/lapack/explore-html/index.html
  !! (Just search the webpage for dgesv or sgesv and you will find a very clear doc)
  !!
  !! Copies x into y
  !!
  !! Interface:
  !!   N    [in]  -> number of elements in vectors
  !!   DX   [in]  -> Source vector. Size > ( 1 + ( N - 1 )*abs( INCX ) )
  !!   INCX [in]  -> Storage spacing between elements of DX
  !!   DY   [out] -> Target vector. Size > ( 1 + ( N - 1 )*abs( INCY ) )
  !!   INCY [in]  -> Storage spacing between elements of DY
  !!
  !! NOTE: It should be assumed that effective size of DX and DY, taking increments under account,
  !!       has to match.
  !!
  interface blas_copy
    !!
    !! Double precision. 64-bit real
    !!
    subroutine dcopy(N, DX, INCX, DY, INCY)
      use iso_fortran_env, only : real64, int32
      implicit none
      integer(int32), intent(in)              :: N
      real(real64), dimension(:), intent(in)  :: DX
      integer(int32), intent(in)              :: INCX
      real(real64), dimension(:), intent(out) :: DY
      integer(int32), intent(in)              :: INCY
    end subroutine dcopy

    !!
    !! Single precision. 32-bit real
    !!
    subroutine scopy(N, SX, INCX, SY, INCY)
      use iso_fortran_env, only : real32, int32
      implicit none
      integer(int32), intent(in)              :: N
      real(real32), dimension(:), intent(in)  :: SX
      integer(int32), intent(in)              :: INCX
      real(real32), dimension(:), intent(out) :: SY
      integer(int32), intent(in)              :: INCY
    end subroutine scopy

  end interface blas_copy


!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Module VARIABLES
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !! This variable must be private to each OpenMP thread
  real(defReal),dimension(:),allocatable,target :: workspace

contains

  !!
  !! Solves linear system of equations of the form Ax=b
  !!
  !! A - any NxN square real matrix
  !! b - real vector of RHS of size N
  !! x - resul vector of size N
  !!
  !! Gives fatalError if input is invalid or solution A is singular
  !!
  subroutine solve(A, x, b)
    real(defReal),dimension(:,:),intent(in) :: A
    real(defReal),dimension(:), intent(out) :: x
    real(defReal),dimension(:), intent(in)  :: b
    real(defReal),dimension(:,:),pointer    :: A_t, B_t
    integer(shortInt),dimension(size(x))    :: pivot
    integer(shortInt)                       :: N, mem, info
    character(100),parameter :: Here='solve ( linearAlgebra_func.f90)'

    ! Verify size of the inputs
    N = size(A,1)
    if(size(b) /= N) then
      call fatalError(Here,'Invalid size of RHS vector b. It is not size N')

    else if(size(x) /=N) then
      call fatalError(Here,'Invallid size of result vector x. It is not size N')

    else if ( any(shape(A) /= N)) then
      call fatalError(Here,'Invalid shape of array A. Is not NxN')

    end if

    ! Calculate memory required and ensure that memory is avalible
    mem = N*N + N
    call getMem(mem)

    ! Associate workspace memory with different variables
    ! Use pointers to change ranks
    A_t(1:N,1:N) => workspace(1:N*N)
    B_t(1:N,1:1) => workspace(N*N+1 : N*N + N)

    ! Copy input
    A_t      = A
    B_t(:,1) = b

    ! Perform calculation
    call lapack_gesv(N, 1, A_t, N, pivot, B_t, N, info)

    if( info < 0) then
      call fatalError(Here,'LINPACK procedure failed with error: '//numToChar(info))

    else if(info > 0) then
      call fatalError(Here,'LINPACK procedure failed. Matrix A is singular.')
    end if

    ! Copy the results out
    x = B_t(:,1)

  end subroutine solve

  !!
  !! Solves a generalised adjoint problem for a GPT response
  !! NOTE: We use <.,.> to denote inner product and ` to denote Adjoint
  !!
  !!  x = M`*x + s
  !!  subject to: <f,x> = 0
  !!
  !! Where M` is adjoint of M and M is singular.
  !! Assumes that <flux,s> = 0 where flux is principal eigenvector of M*flux = 0, which is
  !! a condition for existance of solution x due to Fredholm alternative theorem.
  !! Some for of orthogonalisation f is required for interation to converge.
  !!
  !! A -> Adjoint matrix M`
  !! x -> Solution vector
  !! s -> Generalised adjoint source
  !! f -> Orthogonalisation vector
  !!
  !! Matrix must be NxN and all vectors must be of size N
  !! Gives fatalError if input is invalid or solution fails to convergae
  !!
  !! Maximum number iterations is hardcoded at 1024.
  !! Convergance criterion is <x,x>  smaller then 1.0E-10
  !! Convergance is checkes every 8 iterations.
  !!
  !!
  !! NOTE: For some reasons beyond human comprehension BLAS subroutines do not work in when userd in
  !!       this function. They have been rplaced with Fortran intrinsics, that seem much more
  !!       robust. TODO: Investigate why BLAS fails and hopefully repair it!
  !!
  subroutine solveAdjointProblem(A,x,s,f)
    real(defReal),dimension(:,:), intent(in) :: A
    real(defReal),dimension(:), intent(out)  :: x
    real(defReal),dimension(:), intent(in)   :: s
    real(defReal),dimension(:), intent(in)   :: f
    real(defReal),dimension(:,:),pointer     :: A_t
    real(defReal),dimension(:),pointer       :: x_t, s_t, f_t, w_t
    integer(shortInt)                        :: N, mem, maxIter, i
    real(defReal)                            :: dot
    logical(defBool)                         :: converged
    character(100), parameter :: Here ='solveAdjointProblem ( linearAlgebra_func.f90)'

    ! Verify inputs
    N = size(A,1)
    if(size(s) /= N) then
      call fatalError(Here,'Invalid size of RHS vector s. It is not size N')

    else if(size(x) /=N) then
      call fatalError(Here,'Invalid size of result vector x. It is not size N')

    else if(size(f) /=N) then
      call fatalError(Here,'Invalid size of orthogonalisation vector f. Is is not size N')

    else if ( any(shape(A) /= N)) then
      call fatalError(Here,'Invalid shape of array A. Is not NxN')

    end if

    ! Calculate required work memory and make sure it is available
    mem = N * N + 4 * N
    call getMem(mem)

    ! Associate workspace memory with diffrent wariables
    A_t(1:N,1:N) => workspace(1        : N*N)
    x_t(1:N)     => workspace(N*N+1    : N*N + N)
    w_t(1:N)     => workspace(N*N+N+1  : N*N + 2*N)
    f_t(1:N)     => workspace(N*N+2*N+1: N*N + 3*N)
    s_t(1:N)     => workspace(N*N+3*N+1: N*N + 4*N)

    ! Copy input
    A_t = A
    w_t = s
    s_t = s
    f_t = f
    x_t = ONE

    ! Perform calculation
    maxIter = 1024
    converged = .false.

    do i = 1,maxIter
      ! Calculate w_t = ONE * A_t *x_t + w_t
      w_t = matmul(A_t, x_t) + w_t

      ! Calculate w_t = w_t - <f_t,w_t> * f_t
      dot = dot_product(f_t, w_t) !blas_dot(N, f_t, 1, w_t, 1)
      w_t = w_t - dot * f_t

      ! Normalise solution
      w_t = w_t / norm2(w_t)

      ! Convergance check
      if(iand(i,8) == 0) then
        dot = dot_product(w_t, x_t)
        if (dot < 1.0E-10_defReal) then
          converged = .true.
          exit
        end if
      end if

      ! Move w_t to x_t and s_t to w_t
      x_t = w_t
      w_t = s_t
    end do

    ! Check for failed iteration
    if(.not.converged) then
      call fatalError(Here,'Calculation failed to convarged in: ' // numToChar(maxIter) //' steps. &
      & Final residual is: ' // numToChar(dot))
    end if

    ! Copy out result
    x = w_t

  end subroutine solveAdjointProblem

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
  !! Gives fatalError if input is invalid
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
