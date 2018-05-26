module isotropicAngle_class

  use numPrecision
  use RNG_class,          only : RNG
  use angleLawENDF_inter, only : angleLawENDF
  use isotropicMu_class,  only : isotropicMu


  implicit none
  private

  interface isotropicAngle
    module procedure new_isotropicAngle
  end interface

  !!
  !! Class with mu isotropic at all collisions energies
  !!
  type, public,extends(angleLawENDF) :: isotropicAngle
    private
    type(isotropicMu)   :: muPdf
  contains
    procedure :: sample
    procedure :: probabilityOf
  end type isotropicAngle

contains

  !!
  !! Given collison energy and random number generator sample mu
  !!
  function sample(self,E,rand) result (mu)
    class(isotropicAngle), intent(in) :: self
    real(defReal), intent(in)         :: E
    class(RNG), intent(inout)         :: rand
    real(defReal)                     :: mu

    mu = self % muPdf % sample(rand)

  end function sample

  !!
  !! Return probability density of mu at collision energy E
  !!
  function probabilityOf(self,mu,E) result (prob)
    class(isotropicAngle), intent(in) :: self
    real(defReal), intent(in)         :: E, mu
    real(defReal)                     :: prob

    prob = self % muPdf % probabilityOf(mu)

  end function probabilityOf

  !!
  !! Constructor of isotropicAngle
  !!
  function new_isotropicAngle() result(new)
    type(isotropicAngle) :: new
    ! Nothing to be done
  end function new_isotropicAngle

end module isotropicAngle_class
