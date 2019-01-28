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

  !!


  !!
  !! scoreMemory is a class that stores space for scores from tallies
  !! it is separate from tallyClerks and individual responses to allow:
  !!   -> Easy writing and (later) reading from file for archivisation of results
  !!   -> Easy possibility of extention to tally higher moments of result
  !!   -> Possibility of extention to tally covariance of selected tally bins
  !!   -> Easy copying and recombination of results for OpenMP shared memory parallelism
  !!   -> Easy, output format independent way to perform regression tests
  !!   -> Easy handling of different batch sizes
  !!
  !! For every bin index there are three position, BIN, CSUM, CSUM2. All are initialised to 0.
  !! Function score accumulates value on BIN under given index
  !! Function accumulate accumulates value on CSUM and value^2 on CSUM2 under given index
  !! Subroutine closeCycle, increments cycle counters and possibly closes batch by multiplying all
  !!   scores by a normalisation factor and moving them on the cumulative sum.
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
  !!        bins(binIndex,binType) binType is BIN/CSUM/CSUM2
  !! NOTE2: If batch size is not a denominator of cycles scored results accumulated
  !!        in extra cycles are discarded in current implementation
  !!
  type, public :: scoreMemory
      private
      real(defReal),dimension(:,:),allocatable :: bins          !! Space for binning (2nd dim size is always 3!)
      integer(longInt)                         :: N = 0         !! Size of memory (number of bins)
      integer(shortInt)                        :: id            !! Id of the tally
      integer(shortInt)                        :: batchN = 0    !! Number of Batches
      integer(shortInt)                        :: cycles = 0    !! Cycles counter
      integer(shortInt)                        :: batchSize = 1 !! Batch interval sieze (in cycles)
  contains
    ! Interface procedures
    procedure :: init
    procedure :: kill
    generic   :: score      => score_defReal, score_shortInt, score_longInt
    generic   :: accumulate => accumulate_defReal, accumulate_shortInt, accumulate_longInt
    generic   :: getResult  => getResult_withSTD, getResult_withoutSTD
    procedure :: getScore
    procedure :: closeCycle
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
   self % N = 0
   self % batchN = 0

  end subroutine kill

  !!
  !! Score a result on a given single bin under idx
  !!
  subroutine score_defReal(self, score, idx)
    class(scoreMemory), intent(inout) :: self
    real(defReal), intent(in)         :: score
    integer(longInt), intent(in)     :: idx
    character(100),parameter :: Here = 'score_defReal (scoreMemory_class.f90)'

    ! Verify bounds for the index
    if( idx < 0_longInt .or. idx > self % N) then
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
    integer(longInt), intent(in)     :: idx

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
    integer(longInt), intent(in)     :: idx
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
    integer(longInt), intent(in)     :: idx

    call self % accumulate_defReal(real(score, defReal), idx)

  end subroutine accumulate_longInt

  !!
  !! Close Cylce
  !! Increments cycle counter and detects end-of-batch
  !! When batch finishes it normalises all scores by the factor and moves them to CSUMs
  !!
  subroutine closeCycle(self, normFactor)
    class(scoreMemory), intent(inout) :: self
    real(defReal),intent(in)          :: normFactor

    ! Increment Cycle Counter
    self % cycles = self % cycles + 1

    if(mod(self % cycles, self % batchSize) == 0) then ! Close Batch
      ! Normalise scores
      self % bins(:,BIN) = self % bins(:,BIN) * normFactor

      ! Increment cumulative sums
      self % bins(:,CSUM)  = self % bins(:,CSUM) + self % bins(:,BIN)
      self % bins(:,CSUM2) = self % bins(:,CSUM2) + self % bins(:,BIN) * self % bins(:,BIN)

      ! Zero all score bins
      self % bins(:,BIN) = ZERO

      ! Increment batch counter
      self % batchN = self % batchN + 1

    end if
  end subroutine closeCycle

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
    integer(longInt), intent(in)          :: idx
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

    if(idx < 0_longInt .or. idx > self % N) then
      score = ZERO
    else
      score = self % bins(idx, BIN)
    end if

  end function getScore

end module scoreMemory_class
