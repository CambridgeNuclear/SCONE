module noEnergy_class

  use numPrecision
  use energyLawENDF_inter, only : energyLawENDF
  use RNG_class,           only : RNG

  implicit none
  private

  interface noEnergy
    module procedure new_noEnergy
  end interface

  !!
  !! Null object if for no energy distribution
  !! It is delta function at energy of the indedent particle
  !!
  type, public,extends(energyLawENDF) :: noEnergy
    private
  contains
    procedure :: sample
    procedure :: probabilityOf
  end type noEnergy

contains

  !!
  !! "Sample" outgoing energy.
  !! Returns incedent energy
  !!
  function sample(self,E_in,rand) result (E_out)
    class(noEnergy), intent(in)  :: self
    real(defReal), intent(in)    :: E_in
    class(RNG), intent(inout)    :: rand
    real(defReal)                :: E_out

    E_out = E_in

  end function sample

  !!
  !! Return probability of particle beeing emit with E_out given incedent energy E_in
  !! Returns 1.0 if E_out == E_in (perfectly)
  !!
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

  !!
  !! Constructor
  !!
  function new_noEnergy() result(new)
    type(noEnergy) :: new
    real(defReal) :: dummy
    ! Nothing to be done Avoid compiler warning
    return
    dummy = new % probabilityOf(HALF,ONE)
  end function new_noEnergy


end module noEnergy_class
