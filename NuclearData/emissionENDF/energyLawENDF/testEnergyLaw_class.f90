module testEnergyLaw_class
  use numPrecision
  use RNG_Class,           only : RNG
  use energyLawENDF_inter, only : energyLawENDF

  implicit none
  private

  !!
  !! Energy Law used for unit testing
  !!
  !! Always returns a single value
  !!
  !! Public members:
  !!   E_out -> value returned during sample procedure
  !!
  !! Interface:
  !!   sample        -> returns E_out
  !!   probabilityOf -> returns 1.0 if E_in == E_out (no floating point tolerance!)
  !!
  type, public,extends(energyLawENDF) :: testEnergyLaw
    real(defReal) :: E_out = ZERO
  contains
    ! Interface procedures
    procedure :: sample
    procedure :: probabilityOf
    procedure :: kill
  end type testEnergyLaw

contains

  !!
  !! Sample outgoing energy
  !!
  !! Args:
  !!   E_in [in] -> incident energy [MeV]
  !!   rand [inout] -> random number generator
  !!
  !! Returns:
  !!   Outgoing energy [MeV]
  !!
  function sample(self, E_in, rand) result (E_out)
    class(testEnergyLaw), intent(in) :: self
    real(defReal), intent(in)             :: E_in
    class(RNG), intent(inout)             :: rand
    real(defReal)                         :: E_out

    E_out = self % E_out

  end function sample

  !!
  !! Return probability of transition
  !!
  !! Args:
  !!   E_out [in] -> outgoing energy energy [MeV]
  !!   E_in [in]  -> incident energy [MeV]
  !!
  !! Returns:
  !!   Probability of transition in [0,1]
  !!
  function probabilityOf(self, E_out, E_in) result (prob)
    class(testEnergyLaw), intent(in) :: self
    real(defReal), intent(in)             :: E_out,E_in
    real(defReal)                         :: prob

    if( E_out == self % E_out) then
      prob = ONE
    else
      prob = ZERO
    end if

  end function probabilityOf

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(testEnergyLaw), intent(inout) :: self

    self % E_out = ZERO

  end subroutine kill

end module testEnergyLaw_class
