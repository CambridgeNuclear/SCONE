module isotropicAngle_class

  use numPrecision
  use angleLawENDF_class, only : angleLawENDF
  use isotropicMu_class,  only : isotropicMu
  use RNG_class,          only : RNG

  implicit none
  private

  type, public,extends(angleLawENDF) :: isotropicAngle
    private
    type(isotropicMu)   :: muPdf
  contains
    procedure :: sample
    procedure :: probabilityOf
  end type isotropicAngle

contains

  function sample(self,E,rand) result (mu)
    class(isotropicAngle), intent(in) :: self
    real(defReal), intent(in)         :: E
    class(RNG), intent(inout)         :: rand
    real(defReal)                     :: mu

    mu = self % muPdf % sample(rand)

  end function sample


  function probabilityOf(self,mu,E) result (prob)
    class(isotropicAngle), intent(in) :: self
    real(defReal), intent(in)         :: E, mu
    real(defReal)                     :: prob

    prob = self % muPdf % probabilityOf(mu)

  end function probabilityOf

    
  function new_isotropicAngle()
    type(isotropicAngle),pointer :: new_isotropicAngle

    allocate(new_isotropicAngle)

  end function new_isotropicAngle


end module isotropicAngle_class
