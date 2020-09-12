module testCorrelatedLaw_class
  use numPrecision
  use RNG_class,               only : RNG
  use correlatedLawENDF_inter, only : correlatedLawENDF
  implicit none
  private

  !!
  !! Correlated Law placeholder for testing
  !!
  !! Always returns the same mu and E_out
  !!
  !! Public members:
  !!   mu    -> Returned value of mu
  !!   E_out -> Returned value of E_out
  !!
  !! Interface:
  !!  correlatedLawENDF interface
  !!
  type, public, extends(correlatedLawENDF) :: testCorrelatedLaw
    real(defReal) :: mu    = ZERO
    real(defReal) :: E_out = ZERO
  contains
    ! correlated Law interface
    procedure :: sample
    procedure :: probabilityOf
    procedure :: kill
  end type testCorrelatedLaw

contains

  !!
  !! Sample outgoing angle & energy
  !!
  !! See correlatedLawENDF_inter for detauils
  !!
  subroutine sample(self, mu, E_out, E_in, rand)
    class(testCorrelatedLaw), intent(in) :: self
    real(defReal), intent(out)           :: mu
    real(defReal), intent(out)           :: E_out
    real(defReal), intent(in)            :: E_in
    class(RNG), intent(inout)            :: rand

    mu    = self % mu
    E_out = self % E_out

  end subroutine sample

  !!
  !! Return probability of a given outgoing angle & energy
  !!
  !! See correlatedLawENDF_inter for detauils
  !!
  function probabilityOf(self, mu, E_out, E_in) result(prob)
    class(testCorrelatedLaw), intent(in) :: self
    real(defReal), intent(in)            :: mu
    real(defReal), intent(in)            :: E_out
    real(defReal), intent(in)            :: E_in
    real(defReal)                        :: prob

    if (mu == self % mu .and. E_out == self % E_out ) then
      prob = ONE
    else
      prob = ZERO
    end if

  end function probabilityOf

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(testCorrelatedLaw), intent(inout) :: self

    self % E_out = ZERO
    self % mu = ZERO

  end subroutine kill


end module testCorrelatedLaw_class
