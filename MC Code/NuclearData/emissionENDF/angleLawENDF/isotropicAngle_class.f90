module isotropicAngle_class

  use numPrecision
  use RNG_class,          only : RNG
  use angleLawENDF_inter, only : angleLawENDF
  use aceCard_class,      only : aceCard
  use isotropicMu_class,  only : isotropicMu


  implicit none
  private

  !!
  !! Constructor
  !!
  interface isotropicAngle
    module procedure new_isotropicAngle
  end interface

  !!
  !! Class with mu isotropic at all collisions energies
  !!
  type, public,extends(angleLawENDF) :: isotropicAngle
    private
    type(isotropicMu)   :: muPdf
  contains
    procedure :: init
    procedure :: sample
    procedure :: probabilityOf
    procedure :: kill
  end type isotropicAngle

contains

  !!
  !! Initialise from aceCard and MT number
  !!
  subroutine init(self, ACE, MT)
    class(isotropicAngle), intent(inout):: self
    class(aceCard), intent(inout)       :: ACE
    integer(shortInt), intent(in)       :: MT

    ! Do nothing
    ! No initialisation is needed

  end subroutine init

  !!
  !! Given collison energy and random number generator sample mu
  !!
  function sample(self,E,rand) result (mu)
    class(isotropicAngle), intent(in) :: self
    real(defReal), intent(in)         :: E
    class(RNG), intent(inout)         :: rand
    real(defReal)                     :: mu

    mu = self % muPdf % sample(rand)

  end function sample

  !!
  !! Return probability density of mu at collision energy E
  !!
  function probabilityOf(self,mu,E) result (prob)
    class(isotropicAngle), intent(in) :: self
    real(defReal), intent(in)         :: E, mu
    real(defReal)                     :: prob

    prob = self % muPdf % probabilityOf(mu)

  end function probabilityOf

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(isotropicAngle), intent(inout) :: self

    ! Nothing to do

  end subroutine kill

  !!
  !! Constructor of isotropicAngle
  !!
  function new_isotropicAngle() result(new)
    type(isotropicAngle) :: new
    real(defReal)        :: dummy
    ! Nothing to be done Avoid compiler warning
    return
    dummy = new % probabilityOf(HALF,ONE)
  end function new_isotropicAngle

end module isotropicAngle_class
