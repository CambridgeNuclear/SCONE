module isotropicMu_class

  use numPrecision
  use RNG_class,          only : RNG
  use muEndfPdf_inter,    only : muEndfPdf

  implicit none
  private

  interface isotropicMu
    module procedure new_isotropicMu
  end interface

  !!
  !! Class that contains isotropic secondary angle distribution.
  !! Extends muEndfPdf abstract interface
  !!
  type, public,extends(muEndfPdf) :: isotropicMu
    private
  contains
    procedure :: sample
    procedure :: probabilityOf
  end type isotropicMu

contains

  !!
  !! Samples angle given the random number generator
  !!
  function sample(self,rand) result(mu)
    class(isotropicMu), intent(in)   :: self
    class(RNG), intent(inout)        :: rand
    real(defReal)                    :: mu

    mu = TWO * rand % get() - ONE

  end function sample

  !!
  !! Return probability of mu
  !! Does no check if mu is in <-1;1>
  !!
  function probabilityOf(self,mu) result(probability)
    class(isotropicMu), intent(in)  :: self
    real(defReal), intent(in)       :: mu
    real(defReal)                   :: probability

    probability = HALF

  end function probabilityOf


  !!
  !! Constructor. Needs no arguments
  !!
  function new_isotropicMu() result(new)
    type(isotropicMu):: new
    real(defReal) :: dummy
    ! Nothing to be done Avoid compiler warning
    return
    dummy = new % probabilityOf(HALF)
  end function new_isotropicMu


end module isotropicMu_class
