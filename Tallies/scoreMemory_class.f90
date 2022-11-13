module scoreMemory_class

  use numPrecision
  use genericProcedures, only : fatalError, numToChar
  use openmp_func,       only : ompGetMaxThreads, ompGetThreadNum

  implicit none
  private

  !! Parameters for indexes of per cycle SCORE, Cumulative Sum and Cumulative Sum of squares
  integer(shortInt), parameter :: CSUM  = 1, &
                                  CSUM2 = 2

  !! Size of the 2nd Dimension of bins
  integer(shortInt), parameter :: DIM2 = 2


  !!
  !! scoreMemory is a class that stores space for scores from tallies.
  !! It is separate from tallyClerks and individual responses to allow:
  !!   -> Easy writing and (later) reading from file for archivisation of results
  !!   -> Easy possibility of extention to tally higher moments of result
  !!   -> Possibility of extension to tally covariance of selected tally bins
  !!   -> Easy copying and recombination of results for OpenMP shared memory parallelism
  !!   -> Easy, output format-independent way to perform regression tests
  !!   -> Easy handling of different batch sizes
  !!
  !! For every bin index there are two positions: CSUM, CSUM2. All are initialised to 0.
  !! For scoring, an array is created with dimension (Nbins,nThreads) to mitigate false sharing.
  !! On accumulation, this array adds to the normal bin array.
  !!
  !! Interface:
  !!     init(N,idBS): Initialise with integer size N and integer id. Optional integer Batch Size.
  !!
  !!     kill(): Elemental. Return to uninitialised state.
  !!
  !!     score(score,idx): Score in the bin under idx. FatalError if idx is outside bounds. Score
  !!         is defReal, shortInt or longInt
  !!
  !!     accumulate(score,idx): Accumulate result in cumulative sums in bin under idx. FatalError
  !!         if idx is outside bounds. Score is defReal, shortInt or longInt.
  !!
  !!     getResult(mean, STD, idx, samples): Retrieve mean value and standard deviation of the
  !!         estimate under idx. Use optional samples to specify number of estimates used to
  !!         evaluate mean and STD from default, which is number of batches in score memory.
  !!         STD is optional.
  !!
  !!     getScore(idx): Return current value of score from bin under idx. FatalError if idx is
  !!         outside bounds.
  !!
  !!     closeBin(normFactor,idx): Multiplies score under bin by normFactor and accumulates it in
  !!         cumulative sums. Then sets the bin to zero.
  !!
  !!     closeCycle(normFactor): Multiplies all scores by normFactor and accumulates them in
  !!         cumulative sums. Sets all scors to zero.
  !!
  !!     lastCycle(): Return true if the next call to closeCycle will close a batch.
  !!
  !!     getBatchSize(): Returns number of cycles that constitute a single batch.
  !!
  !! Example use case:
  !!
  !!  do batches=1,20
  !!    do hist=1,10
  !!      call scoreMem % score(hist,1)        ! Score hist (1,10) in bin 1
  !!      call scoreMem % accumulate(hist,2)   ! Accumulate hist in CSUMs of bin 2
  !!    end do
  !!    call scoreMem % closeCycle(ONE)        ! Close batch without normalisation (factor = ONE)
  !!  end do
  !!
  !!  call scoreMem % getResult(mean,STD,1) ! Get result from bin 1 with STD
  !!  call scoreMem % getResult(mean,2,200) ! Get mean from bin 2 assuming 200 samples
  !!
  !! NOTE:  Following indexing is used in bins class member
  !!        bins(binIndex,binType) binType is CSUM/CSUM2
  !! NOTE2: If batch size is not a denominator of cycles scored results accumulated
  !!        in extra cycles are discarded in current implementation
  !!
  type, public :: scoreMemory
      !private
      real(defReal),dimension(:,:),allocatable :: bins          !! Space for storing cumul data (2nd dim size is always 2!)
      real(defReal),dimension(:,:),allocatable :: parallelBins  !! Space for scoring for different threads
      integer(longInt)                         :: N = 0         !! Size of memory (number of bins)
      integer(shortInt)                        :: nThreads = 0  !! Number of threads used for parallelBins
      integer(shortInt)                        :: id            !! Id of the tally
      integer(shortInt)                        :: batchN = 0    !! Number of Batches
      integer(shortInt)                        :: cycles = 0    !! Cycles counter
      integer(shortInt)                        :: batchSize = 1 !! Batch interval size (in cycles)
  contains
    ! Interface procedures
    procedure :: init
    procedure :: kill
    generic   :: score      => score_defReal, score_shortInt, score_longInt
    generic   :: accumulate => accumulate_defReal, accumulate_shortInt, accumulate_longInt
    generic   :: getResult  => getResult_withSTD, getResult_withoutSTD
    procedure :: getScore
    procedure :: closeCycle
    procedure :: closeBin
    procedure :: lastCycle
    procedure :: getBatchSize

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
  !! Optionaly change batchSize from 1 to any +ve number
  !!
  subroutine init(self, N, id, batchSize )
    class(scoreMemory),intent(inout)      :: self
    integer(longInt),intent(in)           :: N
    integer(shortInt),intent(in)          :: id
    integer(shortInt),optional,intent(in) :: batchSize
    character(100), parameter :: Here= 'init (scoreMemory_class.f90)'

    ! Allocate space and zero all bins
    allocate( self % bins(N, DIM2))
    self % bins = ZERO

    self % nThreads = ompGetMaxThreads()

    ! Note the array padding to avoid false sharing
    allocate( self % parallelBins(N + 64, self % nThreads))
    self % parallelBins = ZERO

    ! Save size of memory
    self % N = N

    ! Assign memory id
    self % id = id

    ! Set batchN, cycles and batchSize to default values
    self % batchN    = 0
    self % cycles    = 0
    self % batchSize = 1

    if(present(batchSize)) then
      if(batchSize > 0) then
        self % batchSize = batchSize
      else
        call fatalError(Here,'Batch Size of: '// numToChar(batchSize) //' is invalid')
      end if
    end if

  end subroutine init

  !!
  !! Deallocate memory and return to uninitialised state
  !!
  subroutine kill(self)
   class(scoreMemory), intent(inout) :: self

   if(allocated(self % bins)) deallocate(self % bins)
   if(allocated(self % parallelBins)) deallocate(self % parallelBins)
   self % N = 0
   self % nThreads = 0
   self % batchN = 0

  end subroutine kill

  !!
  !! Score a result on a given single bin under idx
  !!
  subroutine score_defReal(self, score, idx)
    class(scoreMemory), intent(inout) :: self
    real(defReal), intent(in)         :: score
    integer(longInt), intent(in)      :: idx 
    integer(shortInt)                 :: thread_idx
    character(100),parameter :: Here = 'score_defReal (scoreMemory_class.f90)'

    ! Verify bounds for the index
    !if( idx < 0_longInt .or. idx > self % N) then
    !  call fatalError(Here,'Index '//numToChar(idx)//' is outside bounds of &
    !                        & memory with size '//numToChar(self % N))
    !end if

    ! Add the score
    thread_idx = ompGetThreadNum() + 1
    self % parallelBins(idx, thread_idx) = &
            self % parallelBins(idx, thread_idx) + score

  end subroutine score_defReal

  !!
  !! Score a result with shortInt on a given bin under idx
  !!
  subroutine score_shortInt(self, score, idx)
    class(scoreMemory), intent(inout) :: self
    integer(shortInt), intent(in)     :: score
    integer(longInt), intent(in)      :: idx

    call self % score_defReal(real(score, defReal), idx)

  end subroutine score_shortInt

  !!
  !! Score a result with longInt on a given bin under idx
  !!
  subroutine score_longInt(self, score, idx)
    class(scoreMemory), intent(inout) :: self
    integer(longInt), intent(in)      :: score
    integer(longInt), intent(in)     :: idx

    call self % score_defReal(real(score, defReal), idx)

  end subroutine score_longInt

  !!
  !! Increment the result directly on cumulative sums
  !!
  subroutine accumulate_defReal(self, score, idx)
    class(scoreMemory), intent(inout) :: self
    real(defReal), intent(in)         :: score
    integer(longInt), intent(in)      :: idx
    character(100),parameter :: Here = 'accumulate_defReal (scoreMemory_class.f90)'

    ! Verify bounds for the index
    if( idx < 0_longInt .or. idx > self % N) then
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
    integer(longInt), intent(in)     :: idx

    call self % accumulate_defReal(real(score, defReal), idx)

  end subroutine accumulate_shortInt

  !!
  !! Increment the result directly on cumulative sums with longInt score
  !!
  subroutine accumulate_longInt(self, score, idx)
    class(scoreMemory), intent(inout) :: self
    integer(longInt), intent(in)      :: score
    integer(longInt), intent(in)      :: idx

    call self % accumulate_defReal(real(score, defReal), idx)

  end subroutine accumulate_longInt

  !!
  !! Close Cycle
  !! Increments cycle counter and detects end-of-batch
  !! When batch finishes it normalises all scores by the factor and moves them to CSUMs
  !!
  subroutine closeCycle(self, normFactor)
    class(scoreMemory), intent(inout) :: self
    real(defReal),intent(in)          :: normFactor
    integer(longInt)                  :: i
    real(defReal), save               :: res
    !$omp threadprivate(res)

    ! Increment Cycle Counter
    self % cycles = self % cycles + 1

    if(mod(self % cycles, self % batchSize) == 0) then ! Close Batch
      
      !$omp parallel do
      do i = 1, self % N
        
        ! Normalise scores
        self % parallelBins(i,:) = self % parallelBins(i,:) * normFactor
        res = sum(self % parallelBins(i,:))
        
        ! Zero all score bins
        self % parallelBins(i,:) = ZERO
       
        ! Increment cumulative sums 
        self % bins(i,CSUM)  = self % bins(i,CSUM) + res
        self % bins(i,CSUM2) = self % bins(i,CSUM2) + res * res

      end do
      !$omp end parallel do

      ! Increment batch counter
      self % batchN = self % batchN + 1

    end if

  end subroutine closeCycle

  !!
  !! Close Cycle
  !! Multiplies score in bin under idx by normFactor, accumulates it and sets it to zero
  !!
  subroutine closeBin(self, normFactor, idx)
    class(scoreMemory), intent(inout) :: self
    real(defReal),intent(in)          :: normFactor
    integer(longInt), intent(in)      :: idx
    real(defReal)                     :: res
    character(100),parameter :: Here = 'closeBin (scoreMemory_class.f90)'

    ! Verify bounds for the index
    if( idx < 0_longInt .or. idx > self % N) then
      call fatalError(Here,'Index '//numToChar(idx)//' is outside bounds of &
                            & memory with size '//numToChar(self % N))
    end if

    ! Normalise score
    self % parallelBins(idx, :) = self % parallelBins(idx, :) * normFactor

    ! Increment cumulative sum
    res = sum(self % parallelBins(idx,:))
    self % bins(idx,CSUM)  = self % bins(idx,CSUM) + res
    self % bins(idx,CSUM2) = self % bins(idx,CSUM2) + res * res

    ! Zero the score
    self % parallelBins(idx,:) = ZERO

  end subroutine closeBin


  !!
  !! Return true if next closeCycle will close a batch
  !!
  function lastCycle(self) result(isIt)
    class(scoreMemory), intent(in) :: self
    logical(defBool)               :: isIt

    isIt =  mod(self % cycles + 1, self % batchSize) == 0

  end function lastCycle

  !!
  !! Return batchSize
  !!
  pure function getBatchSize(self) result(S)
    class(scoreMemory), intent(in) :: self
    integer(shortInt)              :: S

    S = self % batchSize

  end function getBatchSize

  !!
  !! Load mean result and Standard deviation into provided arguments
  !! Load from bin indicated by idx
  !! Returns 0 if index is invalid
  !!
  elemental subroutine getResult_withSTD(self, mean, STD, idx, samples)
    class(scoreMemory), intent(in)         :: self
    real(defReal), intent(out)             :: mean
    real(defReal),intent(out)              :: STD
    integer(longInt), intent(in)           :: idx
    integer(shortInt), intent(in),optional :: samples
    integer(shortInt)                      :: N
    real(defReal)                          :: inv_N, inv_Nm1

    !! Verify index. Return 0 if not present
    if( idx < 0_longInt .or. idx > self % N) then
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
  elemental subroutine getResult_withoutSTD(self, mean, idx, samples)
    class(scoreMemory), intent(in)         :: self
    real(defReal), intent(out)             :: mean
    integer(longInt), intent(in)           :: idx
    integer(shortInt), intent(in),optional :: samples
    integer(shortInt)                      :: N

    !! Verify index. Return 0 if not present
    if( idx < 0_longInt .or. idx > self % N) then
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

  !!
  !! Obtain value of a score in a bin
  !! Return ZERO for invalid bin address (idx)
  !!
  elemental function getScore(self, idx) result (score)
    class(scoreMemory), intent(in) :: self
    integer(longInt), intent(in)   :: idx
    real(defReal)                  :: score

    if(idx <= 0_longInt .or. idx > self % N) then
      score = ZERO
    else
      score = sum(self % parallelBins(idx, :))
    end if

  end function getScore

end module scoreMemory_class
