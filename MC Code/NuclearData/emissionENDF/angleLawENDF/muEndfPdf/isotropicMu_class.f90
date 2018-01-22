module isotropicMu_class

  use numPrecision
  use RNG_class,          only : RNG
  use muEndfPdf_class,    only : muEndfPdf

  implicit none
  private

  interface isotropicMu
    module procedure new_isotropicMu
  end interface

  type, public,extends(muEndfPdf) :: isotropicMu
    !! Class that contains isotropic secondary angle distribution.
    private
  contains
    procedure :: sample
    procedure :: probabilityOf
  end type isotropicMu

contains

  function sample(self,rand) result(mu)
    class(isotropicMu), intent(in)   :: self
    class(RNG), intent(inout)        :: rand
    real(defReal)                    :: mu

    mu = 2.0 * rand % get() - 1.0
  end function sample


  function probabilityOf(self,mu) result(probability)
    class(isotropicMu), intent(in)  :: self
    real(defReal), intent(in)       :: mu
    real(defReal)                   :: probability

    probability = 0.5
  end function probabilityOf


  function new_isotropicMu()
    type(isotropicMu),pointer :: new_isotropicMu

    allocate(new_isotropicMu)
  end function new_isotropicMu


end module isotropicMu_class
