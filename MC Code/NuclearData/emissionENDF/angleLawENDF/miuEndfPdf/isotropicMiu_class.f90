module isotropicMiu_class

  use numPrecision
  use RNG_class,          only : RNG
  use miuEndfPdf_class,   only : miuEndfPdf

  implicit none
  private

  interface isotropicMiu
    module procedure new_isotropicMiu
  end interface

  type, public,extends(miuEndfPdf) :: isotropicMiu
    !! Class that contains isotropic secondary angle distribution.
    private
  contains
    procedure :: sample
    procedure :: probabilityOf
  end type isotropicMiu

contains

  function sample(self,rand) result(miu)
    class(isotropicMiu), intent(in)  :: self
    class(RNG), intent(inout)        :: rand
    real(defReal)                    :: miu

    miu = 2.0 * rand % get() - 1.0
  end function sample


  function probabilityOf(self,miu) result(probability)
    class(isotropicMiu), intent(in) :: self
    real(defReal), intent(in)       :: miu
    real(defReal)                   :: probability

    probability = 0.5
  end function probabilityOf


  function new_isotropicMiu() !result(newisotropicMiu)
    type(isotropicMiu),pointer :: new_isotropicMiu

    allocate(new_isotropicMiu)
  end function new_isotropicMiu


end module isotropicMiu_class
