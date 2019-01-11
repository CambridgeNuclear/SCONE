module scoreMemory_class

  use numPrecision
  use genericProcedures, only : fatalError, numToChar

  implicit none
  private

  !! Parameters for indexes of per cycle SCORE, Cumulative Sum and Cumulative Sum of squares
  integer(shortInt), parameter :: BIN   = 1, &
                                  CSUM  = 2, &
                                  CSUM2 = 3
  !! Size of the 2nd Dimension of bins
  integer(shortInt), parameter :: DIM2 = 3

  !!
  !! scoreMemory is a class that stores space for scores from tallies
  !! it is separate from tallyClerks and individual responses to allow:
  !!   -> Easy writing and (later) reading from file for archivisation of results
  !!   -> Easy possibility of extention to tally higher moments of result
  !!   -> Possibility of extention to tally covariance of selected tally bins
  !!   -> Easy copying and recombination of results for OpenMP shared memory parallelism
  !!   -> Easy, output format independent way to perform regression tests
  !!
  !! For every bin index there are three position, BIN, CSUM, CSUM2. All are initialised to 0.
  !! Function score accumulates value on BIN under given index
  !! Function accumulate accumulates value on CSUM and value^2 on CSUM2 under given index
  !! Function normalise multiplies all BINs under every index by a factor
  !! Function closeBatch divides all values under BIN by a factor and adds the resulting value on
  !!  CSUM and value^2 on CSUM2.
  !!
  !! Example use case:
  !!
  !!  do batches=1,20
  !!    do hist=1,10
  !!      call scoreMem % score(hist,1)        ! Score hist (1,10) in bin 1
  !!      call scoreMem % accumulate(hist,2)   ! Accumulate hist in CSUMs of bin 2
  !!    end do
  !!    call scoreMem % normalise(factor)      ! Multiply total scores added since last closeBatch by factor
  !!    call scoreMem % closeBatch(batchWgt)   ! Close batch with batchWgt
  !!  end do
  !!  call scoreMem % getResult(mean,STD,1) ! Get result from bin 1 with STD
  !!  call scoreMem % getResult(mean,2,200) ! Get mean from bin 2 assuming 200 samples
  !!
  !! NOTE: Following indexing is used in bins variable
  !!       bins(binIndex,binType) binType is BIN/CSUM/CSUM2
  type, public :: scoreMemory
      private
      real(defReal),dimension(:,:),allocatable :: bins       !! Space for binning (2nd dim size is always 3!)
      integer(shortInt)                        :: N = 0      !! Size of memory (number of bins)
      integer(shortInt)                        :: id         !! Id of the tally
      integer(shortInt)                        :: batchN = 0 !! Number of Batches
  contains
    ! Interface procedures
    procedure :: init
    procedure :: kill
    generic   :: score      => score_defReal, score_shortInt, score_longInt
    generic   :: accumulate => accumulate_defReal, accumulate_shortInt, accumulate_longInt
    generic   :: getResult  => getResult_withSTD, getResult_withoutSTD
    procedure :: normalise
    procedure :: closeBatch

    ! Private procedures
    procedure, private :: score_defReal
    procedure, private :: score_shortInt
    procedure, private :: score_longInt
    procedure, private :: accumulate_defReal
    procedure, private :: accumulate_shortInt
    procedure, private :: accumulate_longInt
    procedure, private :: getResult_withSTD
    procedure, private :: getResult_withoutSTD

  end type scoreMemory

contains

  !!
  !! Allocate space for the bins given number of bins N
  !!
  subroutine init(self, N, id )
    class(scoreMemory), intent(inout) :: self
    integer(shortInt), intent(in)     :: N
    integer(shortInt), intent(in)     :: id

    ! Allocate space and zero all bins
    allocate( self % bins(N, DIM2))
    self % bins = ZERO

    ! Save size of memory
    self % N = N

    ! Assign memory id
    self % id = id

  end subroutine init

  !!
  !! Deallocate memory and return to uninitialised state
  !!
  subroutine kill(self)
   class(scoreMemory), intent(inout) :: self

   deallocate(self % bins)
   self % N = 0
   self % batchN = 0

  end subroutine kill


  !!
  !! Score a result on a given single bin under idx
  !!
  subroutine score_defReal(self, score, idx)
    class(scoreMemory), intent(inout) :: self
    real(defReal), intent(in)         :: score
    integer(shortInt), intent(in)     :: idx
    character(100),parameter :: Here = 'score_defReal (scoreMemory_class.f90)'

    ! Verify bounds for the index
    if( idx < 0 .or. idx > self % N) then
      call fatalError(Here,'Index '//numToChar(idx)//' is outside bounds of &
                            & memory with size '//numToChar(self % N))
    end if

    ! Add the score
    self % bins(idx, BIN) = self % bins(idx, BIN) + score

  end subroutine score_defReal

  !!
  !! Score a result with shortInt on a given bin under idx
  !!
  subroutine score_shortInt(self, score, idx)
    class(scoreMemory), intent(inout) :: self
    integer(shortInt), intent(in)     :: score
    integer(shortInt), intent(in)     :: idx

    call self % score_defReal(real(score, defReal), idx)

  end subroutine score_shortInt

  !!
  !! Score a result with longInt on a given bin under idx
  !!
  subroutine score_longInt(self, score, idx)
    class(scoreMemory), intent(inout) :: self
    integer(longInt), intent(in)      :: score
    integer(shortInt), intent(in)     :: idx

    call self % score_defReal(real(score, defReal), idx)

  end subroutine score_longInt

  !!
  !! Increment the result directly on cumulative sums
  !!
  subroutine accumulate_defReal(self, score, idx)
    class(scoreMemory), intent(inout) :: self
    real(defReal), intent(in)         :: score
    integer(shortInt), intent(in)     :: idx
    character(100),parameter :: Here = 'accumulate_defReal (scoreMemory_class.f90)'

    ! Verify bounds for the index
    if( idx < 0 .or. idx > self % N) then
      call fatalError(Here,'Index '//numToChar(idx)//' is outside bounds of &
                            & memory with size '//numToChar(self % N))
    end if

    ! Add the score
    self % bins(idx, CSUM)  = self % bins(idx, CSUM)  + score
    self % bins(idx, CSUM2) = self % bins(idx, CSUM2) + score * score

  end subroutine accumulate_defReal

  !!
  !! Increment the result directly on cumulative sums with shortInt score
  !!
  subroutine accumulate_shortInt(self, score, idx)
    class(scoreMemory), intent(inout) :: self
    integer(shortInt), intent(in)     :: score
    integer(shortInt), intent(in)     :: idx

    call self % accumulate_defReal(real(score, defReal), idx)

  end subroutine accumulate_shortInt

  !!
  !! Increment the result directly on cumulative sums with longInt score
  !!
  subroutine accumulate_longInt(self, score, idx)
    class(scoreMemory), intent(inout) :: self
    integer(longInt), intent(in)      :: score
    integer(shortInt), intent(in)     :: idx

    call self % accumulate_defReal(real(score, defReal), idx)

  end subroutine accumulate_longInt

  !!
  !! Multiply all scores by a given factor
  !! Does not check for degeneate factor (-ve or 0)!
  !!
  subroutine normalise(self, factor)
    class(scoreMemory), intent(inout) :: self
    real(defReal), intent(in)         :: factor

    ! Multiply all scores
    self % bins(:,BIN) = self % bins(:,BIN) * factor

  end subroutine normalise

  !!
  !! Close batch
  !! Move all scores to the cumulative sums
  !! Scores are divided by a provided divisor
  !! Score bins are zeroed
  !!
  subroutine closeBatch(self, divisor)
    class(scoreMemory), intent(inout) :: self
    real(defReal), intent(in)         :: divisor
    real(defReal)                     :: inv_divisor
    character(100),parameter :: Here = 'closeBatch (scoreMemory.f90)'

    ! Check input
    if( divisor == ZERO) call fatalError(Here,'Divisor is equal to 0!')
    inv_divisor = ONE/divisor

    ! Increment cumulative sums
    self % bins(:,CSUM)  = self % bins(:,CSUM) + self % bins(:,BIN) * inv_divisor
    self % bins(:,CSUM2) = self % bins(:,CSUM2) + self % bins(:,BIN) * self % bins(:,BIN) * inv_divisor * inv_divisor

    ! Zero all score bins
    self % bins(:,BIN) = ZERO

    ! Increment batch counter
    self % batchN = self % batchN + 1

  end subroutine closeBatch

  !!
  !! Load mean result and Standard deviation into provided arguments
  !! Load from bin indicated by idx
  !! Returns 0 if index is invalid
  !!
  subroutine getResult_withSTD(self, mean, STD, idx, samples)
    class(scoreMemory), intent(in)         :: self
    real(defReal), intent(out)             :: mean
    real(defReal),intent(out)              :: STD
    integer(shortInt), intent(in)          :: idx
    integer(shortInt), intent(in),optional :: samples
    integer(shortInt)                      :: N
    real(defReal)                          :: inv_N, inv_Nm1

    !! Verify index. Return 0 if not present
    if( idx < 0 .or. idx > self % N) then
      mean = ZERO
      STD = ZERO
      return
    end if

    ! Check if # of samples is provided
    if( present(samples)) then
      N = samples
    else
      N = self % batchN
    end if

    ! Calculate mean
    mean = self % bins(idx, CSUM) / N

    ! Calculate STD
    inv_N   = ONE / N
    if( N /= 1) then
      inv_Nm1 = ONE / (N - 1)
    else
      inv_Nm1 = ONE
    end if
    STD = self % bins(idx, CSUM2) *inv_N * inv_Nm1 - mean * mean * inv_Nm1
    STD = sqrt(STD)

  end subroutine getResult_withSTD

  !!
  !! Load mean result provided argument
  !! Load from bin indicated by idx
  !! Returns 0 if index is invalid
  !!
  subroutine getResult_withoutSTD(self, mean, idx, samples)
    class(scoreMemory), intent(in)         :: self
    real(defReal), intent(out)             :: mean
    integer(shortInt), intent(in)          :: idx
    integer(shortInt), intent(in),optional :: samples
    integer(shortInt)                      :: N
    real(defReal)                          :: inv_N, inv_Nm1

    !! Verify index. Return 0 if not present
    if( idx < 0 .or. idx > self % N) then
      mean = ZERO
      return
    end if

    ! Check if # of samples is provided
    if( present(samples)) then
      N = samples
    else
      N = self % batchN
    end if

    ! Calculate mean
    mean = self % bins(idx, CSUM) / N

  end subroutine getResult_withoutSTD


end module scoreMemory_class