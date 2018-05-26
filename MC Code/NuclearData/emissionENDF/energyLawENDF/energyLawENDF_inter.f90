module energyLawEndf_inter

  use numPrecision
  use RNG_class, only : RNG

  implicit none
  private

  !!
  !! Abstract interface for diffrent energy distributions
  !!
  type,abstract, public :: energyLawENDF
      private
    contains
      procedure(sample),deferred        :: sample
      procedure(probabilityOf),deferred :: probabilityOf
  end type energyLawENDF


  abstract interface

    !!
    !! Sample outgoing energy given random number generator and incedent energy
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
    function probabilityOf(self,E_out,E_in) result (prob)
      import :: energyLawEndf,&
                defReal
      class(energyLawEndf), intent(in) :: self
      real(defReal), intent(in)        :: E_out,E_in
      real(defReal)                    :: prob
    end function

  end interface
end module energyLawEndf_inter
