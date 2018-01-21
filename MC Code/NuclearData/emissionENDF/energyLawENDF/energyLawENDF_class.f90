module energyLawEndf_class

  use numPrecision
  use RNG_class, only : RNG

  implicit none
  private

  type,abstract, public :: energyLawENDF
      private
    contains
      procedure(sample),deferred        :: sample
      procedure(probabilityOf),deferred :: probabilityOf
  end type energyLawENDF


  abstract interface

    function sample(self,E_in,rand) result (E_out)
      import :: energyLawEndf,&
                defReal,      &
                RNG
      class(energyLawEndf), intent(in) :: self
      real(defReal), intent(in)        :: E_in
      class(RNG), intent(inout)        :: rand
      real(defReal)                    :: E_out
    end function


    function probabilityOf(self,E_out,E_in) result (prob)
      import :: energyLawEndf,&
                defReal
      class(energyLawEndf), intent(in) :: self
      real(defReal), intent(in)        :: E_out,E_in
      real(defReal)                    :: prob
    end function


  end interface


contains


    
end module energyLawEndf_class
