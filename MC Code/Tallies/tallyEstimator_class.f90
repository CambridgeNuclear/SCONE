module tallyEstimator_class

  use numPrecision

  implicit none
  private

  !!
  !! Class to contain a score accumulated during a cycle or history or equvalent
  !! Internal counter is initialised to 0.0 on creation
  !! INTERFACE:
  !!   add(score) -> adds argument score to internal counter score is defReal,sgortInt or longInt
  !!   get()      -> function, returns value(defReal) of internal counter
  !!   reset()    -> resets internal counter to 0.0
  !! NOTES:
  !!   This class exists to hide thread safe addition to a tally score from a user.
  !!   IT IS NOT THREAD SAFE YET.
  !!
  type, public :: tallyScore
    private
    real(defReal)  :: score = ZERO
  contains
    generic   :: add         => add_tallyScore_defReal, &
                                add_tallyScore_shortInt, &
                                add_tallyScore_longInt
    procedure :: get         => get_tallyScore
    procedure :: reset       => reset_tallyScore

    procedure, private :: add_tallyScore_defReal
    procedure, private :: add_tallyScore_shortInt
    procedure, private :: add_tallyScore_longInt

  end type tallyScore

  !!
  !! Class to store an estimate of a quantity provided over many cycles, histories or batches
  !! Is composed and dependent on tallyScore
  !! INTERFACE:
  !!   addEstimate(est)      -> adds extra estimate to internal storage, estimate(est) can be
  !!                            defReal, shortInt or longInt.
  !!   getEstimate(est,STD,N)-> calculates mean estimate and standard deviation of an estimate.
  !!                            est & STD are defReal. N is shortInt. There is no protection
  !!                            against N = 0. Compiler default rules for division by 0 are
  !!                            followed. For N=1 pseudo-variance estimate is returned with N-1
  !!                            term dropped. This is for estetic reasons. This estimate is
  !!                            unreliable.
  !!   reset()               -> resets internal counters to 0.0
  !! NOTES:
  !!   This class exists to wrap cumulative sum and of estimate and estimate^2.
  !!   It is assumed it will be called beween parallel cycles. IS NOT THREAD SAFE
  !!
  type, public :: tallyCounter
    private
    type(tallyScore) :: cSum
    type(tallyScore) :: cSum2

  contains
    generic   :: addEstimate   => addEstimate_tallyCounter_defReal, &
                                  addEstimate_tallyCounter_shortInt, &
                                  addEstimate_tallyCounter_longInt
    generic   :: getEstimate   => getEstimate_tallyCounter_withSTD, &
                                  getEstimate_tallyCounter_withoutSTD
    procedure :: reset         => reset_tallyCounter

    procedure,private :: getEstimate_tallyCounter_withSTD
    procedure,private :: getEstimate_tallyCounter_withoutSTD

    procedure,private :: addEstimate_tallyCounter_defReal
    procedure,private :: addEstimate_tallyCounter_shortInt
    procedure,private :: addEstimate_tallyCounter_longInt

  end type tallyCounter

  !!
  !! Class that combine tallyScore and tally Counter into one entity.
  !! It is used with estimates based on a single integration (i.e. reaction rates)
  !! INTERFACE:
  !!   add(score)   -> adds argument score to internal tallyScore is defReal,sgortInt or longInt
  !!   get()        -> function, returns value(defReal) in internal tallyScore
  !!   closeBatch(N)-> closes batch. Normalises value in internal tallyScore by N.
  !!                   N can be defReal, shortInt and longInt. If no scores were accumulated
  !!                   during a batch (exactly 0.0 value in tallyScore), batch-wise estimate is 0.0
  !!                   even if N=0. Otherwise Compiler default rules for division by 0 are followed.
  !!   getEstimate(est,STD,Nb)-> calculates mean estimate and standard deviation of an estimate.
  !!                             Assumes that Nb batches vere closed before.
  !!                             est & STD are defReal. Nb is shortInt. There is no protection
  !!                             against Nb = 0. Compiler default rules for division by 0 are
  !!                             followed. For Nb=1 pseudo-variance estimate is returned with Nb-1
  !!                             term dropped. This is for estetic reasons. This estimate is
  !!                             unreliable.
  !!  NOTES:
  !!    Class does not store how many batches were closed to save memory and avoid redundunt
  !!    additions in each class instance.
  !!
  type, public :: tallyEstimator
    private
    type(tallyScore)   :: score
    type(tallyCounter) :: estimate
  contains
    generic   :: add         => add_tallyEstimator_defReal, &
                                add_tallyEstimator_shortInt, &
                                add_tallyEstimator_longInt
    procedure :: get         => get_tallyEstimator
    generic   :: closeBatch  => closeBatch_tallyEstimator_defReal, &
                                closeBatch_tallyEstimator_shortInt, &
                                closeBatch_tallyEstimator_longInt
    generic   :: getEstimate => getEstimate_tallyEstimator_withSTD, &
                                getEstimate_tallyEstimator_withoutSTD
    procedure :: reset       => reset_tallyEstimator

    procedure, private :: add_tallyEstimator_defReal
    procedure, private :: add_tallyEstimator_shortInt
    procedure, private :: add_tallyEstimator_longInt

    procedure, private :: closeBatch_tallyEstimator_defReal
    procedure, private :: closeBatch_tallyEstimator_shortInt
    procedure, private :: closeBatch_tallyEstimator_longInt

    procedure, private :: getEstimate_tallyEstimator_withSTD
    procedure, private :: getEstimate_tallyEstimator_withoutSTD

  end type tallyEstimator
contains

!**************************************************************************************************!
!* tallyScore class procedures                                                                    *!
!**************************************************************************************************!
  !!
  !! Add score(defReal) to a score counter
  !!
  elemental subroutine add_tallyScore_defReal(self,score)
    class(tallyScore), intent(inout) :: self
    real(defReal), intent(in)        :: score

    self % score = self % score + score

  end subroutine add_tallyScore_defReal

  !!
  !! Add score(shortInt) to a score counter
  !!
  elemental subroutine add_tallyScore_shortInt(self,score)
    class(tallyScore), intent(inout) :: self
    integer(shortInt), intent(in)    :: score
    real(defReal)                    :: score_tmp

    ! Convert to real
    score_tmp = real(score,defReal)

    ! Call procedure for real
    call self % add_tallyScore_defReal(score_tmp)

  end subroutine add_tallyScore_shortInt

  !!
  !! Add score(longInt) to a score counter
  !!
  elemental subroutine add_tallyScore_longInt(self,score)
    class(tallyScore), intent(inout) :: self
    integer(longInt), intent(in)     :: score
    real(defReal)                    :: score_tmp

    ! Convert to real
    score_tmp = real(score,defReal)

    ! Call procedure for real
    call self % add_tallyScore_defReal(score_tmp)

  end subroutine add_tallyScore_longInt

  !!
  !! Returns accumulated score value
  !!
  elemental function get_tallyScore(self) result(score)
    class(tallyScore), intent(in) :: self
    real(defReal)                 :: score

    score = self % score

  end function get_tallyScore

  !!
  !! Resets score counter
  !!
  elemental subroutine reset_tallyScore(self)
    class(tallyScore), intent(inout) :: self

    self % score = ZERO

  end subroutine reset_tallyScore

!**************************************************************************************************!
!* tallyCounter class procedures                                                                    *!
!**************************************************************************************************!
  !!
  !! Add a result estimate(defReal) to a tallyCounter
  !!
  elemental subroutine addEstimate_tallyCounter_defReal(self,score)
    class(tallyCounter), intent(inout) :: self
    real(defReal),intent(in)           :: score

    call self % cSum  % add(score)
    call self % cSum2 % add(score * score )

  end subroutine addEstimate_tallyCounter_defReal

  !!
  !! Add a result estimate(shortInt) to a tallyCounter
  !!
  elemental subroutine addEstimate_tallyCounter_shortInt(self,score)
    class(tallyCounter), intent(inout) :: self
    integer(shortInt), intent(in)      :: score
    real(defReal)                      :: score_tmp

    ! Convert and call procedure for defReal
    score_tmp = real(score,defReal)
    call self % addEstimate(score_tmp)

  end subroutine addEstimate_tallyCounter_shortInt

  !!
  !! Add a result estimate(longInt) to a tallyCounter
  !!
  elemental subroutine addEstimate_tallyCounter_longInt(self,score)
    class(tallyCounter), intent(inout) :: self
    integer(longInt), intent(in)       :: score
    real(defReal)                      :: score_tmp

    ! Convert and call procedure for defReal
    score_tmp = real(score,defReal)
    call self % addEstimate(score_tmp)

  end subroutine addEstimate_tallyCounter_longInt

  !!
  !! Obtain estimate from a tally Counter. Returns standard deviation
  !! Normalise assuming N (shortInt) samples were collected
  !! There is no protection against -ve and N = 0 normal float arthmetic rules are followed
  !!
  elemental subroutine getEstimate_tallyCounter_withSTD(self,est,STD,N)
    class(tallyCounter), intent(in) :: self
    real(defReal), intent(out)      :: est
    real(defReal), intent(out)      :: STD
    integer(shortInt), intent(in)   :: N
    real(defReal)                   :: cSum2
    integer(shortInt)               :: Nm1

    cSum2  = self % cSum2  % get()

    ! Call getEstimate without STD to calculate mean
    call self % getEstimate(est,N)

    ! Protect against STD = NaN if N == 1. Return variance equal to estimate
    ! This mostly done for aesthetics. Single sample estimate is unrealible
    if ( N == 1) then
      STD = est

    else ! when N /= 1
      ! Precalculate denominator (N-1)
      Nm1 = N - 1

      ! Calculate the Variance of the mean
      STD = cSum2 / Nm1 /N - est * est / Nm1

      ! Calculate Standard Deviation
      STD = sqrt(STD)

    end if

  end subroutine getEstimate_tallyCounter_withSTD

  !!
  !! Obtain estimate from a tally Counter. Does not returns Standard Deviation.
  !! Normalise assuming N (shortInt) samples were collected
  !! There is no protection against -ve and N = 0 normal float arthmetic rules are followed
  !!
  elemental subroutine getEstimate_tallyCounter_withoutSTD(self,est,N)
    class(tallyCounter), intent(in) :: self
    real(defReal), intent(out)      :: est
    integer(shortInt), intent(in)   :: N
    real(defReal)                   :: cSum

    cSum  = self % cSum  % get()

    ! Calculate Mean
    est = cSum / N

  end subroutine getEstimate_tallyCounter_withoutSTD

  !!
  !! Reset tallyCounter
  !!
  elemental subroutine reset_tallyCounter(self)
    class(tallyCounter), intent(inout) :: self

    call self % cSum  % reset()
    call self % cSum2 % reset()

  end subroutine reset_tallyCounter

!**************************************************************************************************!
!* tallyEstimate class procedures                                                                    *!
!**************************************************************************************************!

  !!
  !! Add new score(defReal) to a tallyScore in the estimator
  !!
  elemental subroutine add_tallyEstimator_defReal(self,score)
    class(tallyEstimator), intent(inout) :: self
    real(defReal), intent(in)            :: score

    call self % score % add(score)

  end subroutine add_tallyEstimator_defReal

  !!
  !! Add new score(defReal) to a tallyScore in the estimator
  !!
  elemental subroutine add_tallyEstimator_shortInt(self,score)
    class(tallyEstimator), intent(inout) :: self
    integer(shortInt), intent(in)        :: score
    real(defReal)                        :: score_tmp

    ! Convert and call procedure for defReal
    score_tmp = real(score,defReal)
    call self % score % add(score)

  end subroutine add_tallyEstimator_shortInt

  !!
  !! Add new score(defReal) to a tallyScore in the estimator
  !!
  elemental subroutine add_tallyEstimator_longInt(self,score)
    class(tallyEstimator), intent(inout) :: self
    integer(longInt), intent(in)        :: score
    real(defReal)                        :: score_tmp

    ! Convert and call procedure for defReal
    score_tmp = real(score,defReal)
    call self % score % add(score)

  end subroutine add_tallyEstimator_longInt

  !!
  !! Return accumulated score
  !!
  elemental function get_tallyEstimator(self) result (score)
    class(tallyEstimator), intent(in) :: self
    real(defReal)                     :: score

    score = self % score % get()

  end function get_tallyEstimator

  !!
  !! Closes Batch
  !! Normalise score assuming that total weight of the samples was N(defReal)
  !! There is no protection against N == 0
  !! normal float arthmetic rules are followed
  !! Special case: if score == 0.0 (exactly in floating points)
  !!   It corresponds (most likley) to no samples beeing collected
  !!   Then batch-wise mean is estimated to be 0.0
  !!
  elemental subroutine closeBatch_tallyEstimator_defReal(self,N)
    class(tallyEstimator), intent(inout) :: self
    real(defReal), intent(in)            :: N
    real(defReal)   :: score, mean

    ! Extract score and reset tallyScore
    score = self % score % get()
    call self % score % reset()

    ! Detect spetial case when no samples were collected
    if( score == ZERO ) then
      mean = ZERO

    else
      mean = score / N

    end if

    ! Add estimate of the mean to the tallyCounter
    call self % estimate % addEstimate(mean)

  end subroutine

  !!
  !! Normalise score assuming that total weight of the samples was N(shortInt)
  !! Converts N into defReal and calls procedure for defReal
  !!
  elemental subroutine closeBatch_tallyEstimator_shortInt(self,N)
    class(tallyEstimator), intent(inout) :: self
    integer(shortInt), intent(in)        :: N
    real(defReal)                        :: N_tmp

    ! Convert and call procedure for defReal
    N_tmp = real(N,defReal)
    call self % closeBatch(N_tmp)

  end subroutine closeBatch_tallyEstimator_shortInt

  !!
  !! Normalise score assuming that total weight of the samples was N(longInt)
  !! Converts N into defReal and calls procedure for defReal
  !!
  elemental subroutine closeBatch_tallyEstimator_longInt(self,N)
    class(tallyEstimator), intent(inout) :: self
    integer(longInt), intent(in)         :: N
    real(defReal)                        :: N_tmp

    ! Convert and call procedure for defReal
    N_tmp = real(N,defReal)
    call self % closeBatch(N_tmp)

  end subroutine closeBatch_tallyEstimator_longInt

  !!
  !! Obtain estimate from a tally Counter. Returns Standard Deviation.
  !! Normalise assuming N (shortInt) batches were collected (with closeBatch)
  !! There is no protection against N < 0 and N == 0
  !! Normal float arthmetic rules are followed
  !!
  elemental subroutine getEstimate_tallyEstimator_withSTD(self,est,STD,Nb)
    class(tallyEstimator), intent(in)    :: self
    real(defReal), intent(out)           :: est
    real(defReal), intent(out)           :: STD
    integer(shortInt), intent(in)        :: Nb

    call self % estimate % getEstimate(est,STD,Nb)

  end subroutine getEstimate_tallyEstimator_withSTD

  !!
  !! Obtain estimate from a tally Counter. Does not return Standard Deviation.
  !! Normalise assuming N (shortInt) batches were collected (with closeBatch)
  !! There is no protection against N < 0 and N == 0
  !! Normal float arthmetic rules are followed
  !!
  elemental subroutine getEstimate_tallyEstimator_withoutSTD(self,est,Nb)
    class(tallyEstimator), intent(in)    :: self
    real(defReal), intent(out)           :: est
    integer(shortInt), intent(in)        :: Nb

    call self % estimate % getEstimate(est,Nb)

  end subroutine getEstimate_tallyEstimator_withoutSTD

  !!
  !! Resets score and estimate of tallyEstimate
  !!
  elemental subroutine reset_tallyEstimator(self)
    class(tallyEstimator), intent(inout) :: self

    call self % score % reset()
    call self % estimate % reset()

  end subroutine reset_tallyEstimator

end module tallyEstimator_class
