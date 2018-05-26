module angleLawENDF_inter

  use numPrecision
  use RNG_class, only : RNG

  implicit none
  private

  !!
  !! Abstract interface for object containing angle (mu) pdf
  !!
  type,abstract, public :: angleLawENDF
    private
  contains
    procedure(sample),deferred         :: sample
    procedure(probabilityOf),deferred  :: probabilityOf
  end type angleLawENDF


  abstract interface

    !!
    !! Given collison energy and random number generator sample mu
    !!
    function sample(self,E,rand) result (mu)
      import :: angleLawENDF, &
                RNG,          &
                defReal
      class(angleLawENDF), intent(in)  :: self
      real(defReal), intent(in)        :: E
      class(RNG), intent(inout)        :: rand
      real(defReal)                    :: mu
    end function sample

    !!
    !! Return probability density of mu at collision energy E
    !!
    function probabilityOf(self,mu,E) result (prob)
      import :: angleLawENDF, &
                defReal
      class(angleLawENDF), intent(in) :: self
      real(defReal), intent(in)       :: E, mu
      real(defReal)                   :: prob
    end function probabilityOf

  end interface
    
end module angleLawENDF_inter
