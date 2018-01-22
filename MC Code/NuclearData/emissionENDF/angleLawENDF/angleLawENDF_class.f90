module angleLawENDF_class

  use numPrecision
  use RNG_class, only : RNG

  implicit none
  private

  type,abstract, public :: angleLawENDF
    private
  contains
    procedure(sample),deferred         :: sample
    procedure(probabilityOf),deferred  :: probabilityOf
  end type angleLawENDF


  abstract interface

    function sample(self,E,rand) result (mu)
      import :: angleLawENDF, &
                RNG,          &
                defReal
      class(angleLawENDF), intent(in)  :: self
      real(defReal), intent(in)        :: E
      class(RNG), intent(inout)        :: rand
      real(defReal)                    :: mu
    end function sample


    function probabilityOf(self,mu,E) result (prob)
      import :: angleLawENDF, &
                defReal
      class(angleLawENDF), intent(in) :: self
      real(defReal), intent(in)       :: E, mu
      real(defReal)                   :: prob
    end function probabilityOf

  end interface

contains
    
end module angleLawENDF_class
