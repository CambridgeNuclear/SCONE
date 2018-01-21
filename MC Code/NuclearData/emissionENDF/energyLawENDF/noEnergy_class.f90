module noEnergy_class

  use numPrecision
  use energyLawENDF_class, only : energyLawENDF
  use RNG_class,           only : RNG

  implicit none
  private

  interface noEnergy
    module procedure new_noEnergy
  end interface

  type, public,extends(energyLawENDF) :: noEnergy
    private
  contains
    procedure :: sample
    procedure :: probabilityOf
  end type noEnergy

contains

  function sample(self,E_in,rand) result (E_out)
    class(noEnergy), intent(in)  :: self
    real(defReal), intent(in)    :: E_in
    class(RNG), intent(inout)    :: rand
    real(defReal)                :: E_out

    E_out = E_in

  end function sample
    
  function probabilityOf(self,E_out,E_in) result (prob)
    class(noEnergy), intent(in)       :: self
    real(defReal), intent(in)         :: E_out, E_in
    real(defReal)                     :: prob

    if (E_out == E_in) then
      prob = 1.0
    else
      prob = 0.0
    end if
  end function probabilityOf

  function new_noEnergy()
    type(noEnergy),pointer :: new_noEnergy

    allocate(new_noEnergy)

  end function


end module noEnergy_class
