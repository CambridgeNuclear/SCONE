module angleLawENDF_inter

  use numPrecision
  use RNG_class,     only : RNG
  use aceCard_class, only : aceCard

  implicit none
  private

  !!
  !! Abstract interface for object containing angle (mu) pdf
  !!
  type,abstract, public :: angleLawENDF
    private
  contains
    procedure(init),deferred           :: init
    procedure(sample),deferred         :: sample
    procedure(probabilityOf),deferred  :: probabilityOf
    procedure(kill),deferred           :: kill
  end type angleLawENDF


  abstract interface

    !!
    !! Initialise angular law from aceCard and MT number
    !!
    subroutine init(self, ACE, MT)
      import :: angleLawENDF, &
                aceCard, &
                shortInt
      class(angleLawENDF), intent(inout) :: self
      class(aceCard), intent(inout)      :: ACE
      integer(shortInt), intent(in)      :: MT
    end subroutine init

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

    !!
    !! Return to uninitialised state
    !!
    elemental subroutine kill(self)
      import :: angleLawENDF
      class(angleLawENDF), intent(inout) :: self
    end subroutine kill


  end interface
    
end module angleLawENDF_inter
