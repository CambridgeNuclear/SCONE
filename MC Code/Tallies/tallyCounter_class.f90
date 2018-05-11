module tallyCounter_class

  use numPrecision
  use genericProcedures, only : fatalError

  implicit none
  private

  !!
  !! Working horse of tallies
  !! Class that handles accumulation of scores
  !! Supports batch scheme to mitigate correlation between cycles
  !! Following subroutines are avalibe:
  !!   score(value)         ->  adds value to the current score
  !!   getScore(mean,STD,N) ->  caluclates mean and standard deviation of the score assuming N
  !!                            samples were stored. If N == 0 or N < 0 returns an error.
  !!   closeBatch(N)        ->  calculates mean and moves it to batch & batch2 component. Resets
  !!                            score and score2 to 0.0. Assumes N is number of samples stored.
  !!                            If N == 0 or N < 0 returns an error.
  !!   getEstimate(Nb)      ->  calculates mean and standard deviation of all batches assuming Nb
  !!                            is a number of all batches stored. If Nb == 0 or Nb < 0 returns an
  !!                            error.
  !!   reset()              -> returns score, score2, batch, batch2 back to 0.0
  !!
  !!  N or Nb can be defReal, shortInt or longInt
  !!
  type, public :: tallyCounter
    private
    real(defReal) :: score  = 0.0 !! Cumulative sum of scores
    real(defReal) :: score2 = 0.0 !! Cumulative sum of scores squared
    real(defReal) :: batch  = 0.0 !! Cumulative sum of batch-wise scores
    real(defReal) :: batch2 = 0.0 !! Cumulative sum of batch-wise scores squared

  contains
    procedure :: addScore
    generic   :: getScore    => getScore_defReal,    getScore_shortInt,    getScore_longInt
    generic   :: closeBatch  => closeBatch_defReal,  closeBatch_shortInt,  closeBatch_longInt
    generic   :: getEstimate => getEstimate_defReal, getEstimate_shortInt, getEstimate_longInt
    procedure :: reset

    ! Local bindings for generic type-boud procedures
    procedure,private :: getScore_defReal
    procedure,private :: getScore_shortInt
    procedure,private :: getScore_longInt

    procedure,private :: closeBatch_defReal
    procedure,private :: closeBatch_shortInt
    procedure,private :: closeBatch_longInt

    procedure,private :: getEstimate_defReal
    procedure,private :: getEstimate_shortInt
    procedure,private :: getEstimate_longInt
  end type tallyCounter

contains

  !!
  !! Add new score to cumulative sums in a tally counter
  !!
  subroutine addScore(self,s)
    class(tallyCounter), intent(inout) :: self
    real(defReal), intent(in)          :: s

    ! Update internal cumulative sums
    self % score  = self % score + s
    self % score2 = self % score2 + s * s

  end subroutine addScore


  !!
  !! Calculate mean and standard deviation of the score given number of samples(defReal)
  !!
  subroutine getScore_defReal(self,mean,STD,N)
    class(tallyCounter), intent(in)    :: self
    real(defReal), intent(out)         :: mean
    real(defReal), intent(out)         :: STD
    real(defReal), intent(in)          :: N
    character(100),parameter           :: Here = 'getScore_defReal (tallyCounter_class.f90)'

    ! Proctect again N == 0.0 or N < 0.0
    if ( N <= ZERO ) then
      call fatalError(Here,'Attempted calculation of score with not +ve number of samples')
    end if

    ! Calculate mean
    mean = self % score / N

    ! Calculate Sample Variance
    STD  = (self % score2 -  mean*mean) / (N-ONE)

    ! Calculate Variance of the mean
    STD  = STD / N

    ! Calculate Standard Deviation
    STD = sqrt( STD )

  end subroutine getScore_defReal

  !!
  !! Calculate mean and standard deviation of the score given number of samples(shortInt)
  !!
  subroutine getScore_shortInt(self,mean,STD,N)
    class(tallyCounter), intent(in)    :: self
    real(defReal), intent(out)         :: mean
    real(defReal), intent(out)         :: STD
    integer(shortInt),intent(in)       :: N
    real(defReal)                      :: N_real

    ! Convert integer to real and call procedure to for real to avoid code repeat
    N_real = real(N,defReal)
    call self % getScore_defReal(mean,STD,N_real)

  end subroutine getScore_shortInt

  !!
  !! Calculate mean and standard deviation of the score given number of samples(longInt)
  !!
  subroutine getScore_longInt(self,mean,STD,N)
    class(tallyCounter), intent(in)    :: self
    real(defReal), intent(out)         :: mean
    real(defReal), intent(out)         :: STD
    integer(longInt),intent(in)        :: N
    real(defReal)                      :: N_real

    ! Convert integer to real and call procedure to for real to avoid code repeat
    N_real = real(N,defReal)
    call self % getScore_defReal(mean,STD,N_real)

  end subroutine getScore_longInt


  !!
  !! Add current score to batch-wise estimate given number of samples(N) since last batch (defReal)
  !!
  subroutine closeBatch_defReal(self,N)
    class(tallyCounter), intent(inout) :: self
    real(defReal), intent(inout)       :: N
    real(defReal)                      :: mean
    character(100),parameter           :: Here = 'getScore_defReal (tallyCounter_class.f90)'

    ! Proctect again N == 0.0 or N < 0.0
    if ( N <= ZERO ) then
      call fatalError(Here,'Attempted calculation of score with not +ve number of samples')
    end if

    ! Calculate mean
    mean = self % score / N

    ! Add estimate of the mean from this batch to cumulative batch-wise scores
    self % batch  = self % batch  + mean
    self % batch2 = self % batch2 + mean * mean

    ! Reset score counters
    self % score  = ZERO
    self % score2 = ZERO

  end subroutine closeBatch_defReal

  !!
  !! Add current score to batch-wise estimate given number of samples(N) since last batch (defReal)
  !!
  subroutine closeBatch_shortInt(self,N)
    class(tallyCounter), intent(inout) :: self
    integer(shortInt),intent(in)       :: N
    real(defReal)                      :: N_real

    ! Convert integer to real and call procedure to for real to avoid code repeat
    N_real = real(N,defReal)
    call self % closeBatch_defReal(N_real)

  end subroutine closeBatch_shortInt

  !!
  !! Add current score to batch-wise estimate given number of samples(N) since last batch (defReal)
  !!
  subroutine closeBatch_longInt(self,N)
    class(tallyCounter), intent(inout) :: self
    integer(longInt),intent(in)        :: N
    real(defReal)                      :: N_real

    ! Convert integer to real and call procedure to for real to avoid code repeat
    N_real = real(N,defReal)
    call self % closeBatch_defReal(N_real)

  end subroutine closeBatch_longInt

  !!
  !! Calculate mean and standard deviation of the result estimate
  !! given number of batches Nb (defReal)
  !!
  subroutine getEstimate_defReal(self,mean,STD,Nb)
    class(tallyCounter), intent(inout) :: self
    real(defReal), intent(out)         :: mean
    real(defReal), intent(out)         :: STD
    real(defReal), intent(in)          :: Nb
    character(100), parameter          :: Here = 'getEstimate_defReal (tallyCounter_class.f90)'

    ! Proctect again Nb == 0.0 or Nb < 0.0
    if ( Nb <= ZERO ) then
      call fatalError(Here,'Attempted calculation of estimate with not +ve number of samples')
    end if

    ! Calculate mean
    mean = self % batch / Nb

    ! Calculate Sample Variance
    STD  = (self % batch2 -  mean*mean) / (Nb-ONE)

    ! Calculate Variance of the mean
    STD  = STD / Nb

    ! Calculate Standard Deviation
    STD = sqrt( STD )

  end subroutine getEstimate_defReal

  !!
  !! Calculate mean and standard deviation of the result estimate
  !! given number of batches Nb (shortInt)
  !!
  subroutine getEstimate_shortInt(self,mean,STD,Nb)
    class(tallyCounter), intent(inout) :: self
    real(defReal), intent(out)         :: mean
    real(defReal), intent(out)         :: STD
    integer(shortInt),intent(in)       :: Nb
    real(defReal)                      :: Nb_real

    ! Convert integer to real and call procedure to for real to avoid code repeat
    Nb_real = real(Nb,defReal)
    call self % getEstimate_defReal(mean,STD,Nb_real)


  end subroutine getEstimate_shortInt

  !!
  !! Calculate mean and standard deviation of the result estimate
  !! given number of batches Nb (longInt)
  !!
  subroutine getEstimate_longInt(self,mean,STD,Nb)
    class(tallyCounter), intent(inout) :: self
    real(defReal), intent(out)         :: mean
    real(defReal), intent(out)         :: STD
    integer(longInt),intent(in)        :: Nb
    real(defReal)                      :: Nb_real

    ! Convert integer to real and call procedure to for real to avoid code repeat
    Nb_real = real(Nb,defReal)
    call self % getEstimate_defReal(mean,STD,Nb_real)

  end subroutine getEstimate_longInt

  !!
  !! Sets tallyCounter back to initial state (all components are 0.0)
  !!
  subroutine reset(self)
    class(tallyCounter), intent(inout) :: self

    ! Reset all components
    self % score  = ZERO
    self % score2 = ZERO
    self % batch  = ZERO
    self % batch2 = ZERO

  end subroutine reset

end module tallyCounter_class
