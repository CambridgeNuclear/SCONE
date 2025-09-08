module noAngle_class

  use numPrecision
  use RNG_class,          only : RNG
  use aceCard_class,      only : aceCard
  use angleLawENDF_inter, only : angleLawENDF


  implicit none
  private

  interface noAngle
    module procedure new_noAngle
  end interface

  !!
  !! Null object for lack of mu data
  !! Can be used as a placeholder for reactions that are no scattering
  !!
  type, public,extends(angleLawENDF) :: noAngle
    private
  contains
    procedure :: init
    procedure :: sample
    procedure :: probabilityOf
    procedure :: kill
  end type noAngle

contains

  !!
  !! Initialise from aceCard and MT number
  !!
  subroutine init(self, ACE, MT)
    class(noAngle), intent(inout)  :: self
    class(aceCard), intent(inout)  :: ACE
    integer(shortInt), intent(in)  :: MT

    ! Do nothing
    ! No initialisation is needed

  end subroutine init


  !!
  !! Do not change direction
  !! Always return mu = 1.0
  !!
  function sample(self,E,rand) result (mu)
    class(noAngle), intent(in)   :: self
    real(defReal), intent(in)    :: E
    class(RNG), intent(inout)    :: rand
    real(defReal)                :: mu

    mu = 1.0

  end function sample

  !!
  !! Does not change direction
  !! Probability density is a delta function delta(mu-1.0)
  !!
  function probabilityOf(self,mu,E) result (prob)
    class(noAngle), intent(in)  :: self
    real(defReal), intent(in)   :: E, mu
    real(defReal)               :: prob

    if (mu == 1.0_defReal) then
      prob = 1.0
    else
      prob = 0.0
    end if

  end function probabilityOf

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(noAngle), intent(inout) :: self

    ! Nothing to do

  end subroutine kill

  !!
  !! Constructor
  !!
  function new_noAngle() result(new)
    type(noAngle)  :: new
    real(defReal) :: dummy
    ! Nothing to be done Avoid compiler warning
    return
    dummy = new % probabilityOf(HALF,ONE)
  end function new_noAngle


end module noAngle_class
