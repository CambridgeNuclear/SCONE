module energyLawEndf_inter

  use numPrecision
  use RNG_class, only : RNG

  implicit none
  private

  !!
  !! Abstract interface for diffrent energy distributions
  !!
  !! Interface:
  !!   sample        -> Sample outgoing energy
  !!   probabilityOf -> Return propability density at outgoing the energy & angle
  !!   kill          -> Return to uninitialised state
  !!
  type,abstract, public :: energyLawENDF
      private
    contains
      procedure(sample),deferred        :: sample
      procedure(probabilityOf),deferred :: probabilityOf
      procedure(kill),deferred          :: kill
  end type energyLawENDF


  abstract interface

    !!
    !! Sample outgoing energy given random number generator and incedent energy
    !!
    !! Args:
    !!   E_in [in] -> incident energy [MeV]
    !!   rand [inout] -> random number generator
    !!
    !! Returns:
    !!   Outgoing energy [MeV]
    !!
    !! Errors:
    !!   Gives fatalError if sampling fails to sample energy Law
    !!
    function sample(self,E_in,rand) result (E_out)
      import :: energyLawEndf,&
                defReal,      &
                RNG
      class(energyLawEndf), intent(in) :: self
      real(defReal), intent(in)        :: E_in
      class(RNG), intent(inout)        :: rand
      real(defReal)                    :: E_out
    end function

    !!
    !! Give probability of outgoing energy given incedent energy
    !!
    !! Args:
    !!   E_out [in] -> outgoing energy energy [MeV]
    !!   E_in [in] -> incident energy [MeV]
    !!
    !! Returns:
    !!   Probability of transition in [0,1]
    !!
    function probabilityOf(self,E_out,E_in) result (prob)
      import :: energyLawEndf,&
                defReal
      class(energyLawEndf), intent(in) :: self
      real(defReal), intent(in)        :: E_out,E_in
      real(defReal)                    :: prob
    end function

    !!
    !! Return to uninitialised state
    !!
    elemental subroutine kill(self)
      import :: energyLawEndf
      class(energyLawEndf), intent(inout) :: self
    end subroutine kill

  end interface
end module energyLawEndf_inter
