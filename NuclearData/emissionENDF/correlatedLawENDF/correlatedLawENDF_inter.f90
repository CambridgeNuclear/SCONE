module correlatedLawENDF_inter

  use numPrecision
  use RNG_class,    only : RNG

  implicit none
  private

  !!
  !! Abstract interface for objects containing correleated mu energy data
  !!
  type, public,abstract :: correlatedLawENDF
    private
    contains
      procedure(sample),deferred        :: sample
      procedure(probabilityOf),deferred :: probabilityOf
      procedure(kill),deferred          :: kill
  end type correlatedLawENDF

  abstract interface

    !!
    !! Samples mu and E_out givent incident energy E_in and random nummber generator
    !!
    subroutine sample(self,mu,E_out,E_in,rand)
      import :: correlatedLawENDF, &
                defReal, &
                RNG
      class(correlatedLawENDF), intent(in) :: self
      real(defReal), intent(out)           :: mu
      real(defReal), intent(out)           :: E_out
      real(defReal), intent(in)            :: E_in
      class(RNG), intent(inout)            :: rand
    end subroutine

    !!
    !! Returns probability that neutron was emmited at mu & E_out given incident energy E_in
    !!
    function probabilityOf(self,mu,E_out,E_in) result(prob)
      import :: correlatedLawENDF, &
                defReal
      class(correlatedLawENDF), intent(in) :: self
      real(defReal), intent(in)            :: mu
      real(defReal), intent(in)            :: E_out
      real(defReal), intent(in)            :: E_in
      real(defReal)                        :: prob
    end function probabilityOf

    !!
    !! Return to uninitialised state
    !!
    elemental subroutine kill(self)
      import :: correlatedLawENDF
      class(correlatedLawENDF), intent(inout) :: self
    end subroutine kill

  end interface

end module correlatedLawENDF_inter
