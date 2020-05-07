module polynomialrelease_class

  use numPrecision
  use aceCard_class,        only : aceCard
  use releaseLawENDF_inter, only : releaseLawENDF

  implicit none
  private

  interface polynomialRelease
    module procedure new_polynomialRelease
    module procedure new_polynomialRelease_fromACE
  end interface

  !!
  !! Polynomial representation of Nu data
  !! a_0 + a_1 * x + a_2 * x**2 + ... etc.
  !!
  type, public,extends(releaseLawENDF) :: polynomialRelease
      private
      real(defReal),dimension(:),allocatable :: coeffs  !! Polynomial coefficients [a_0,a_1,...]
    contains
      procedure :: init
      procedure :: releaseAt
      procedure :: kill
  end type polynomialRelease

contains

  !!
  !! Initialise
  !!
  subroutine init(self,coeffs)
    class(polynomialRelease), intent(inout)  :: self
    real(defReal),dimension(:),intent(in)    :: coeffs

    if (allocated(self % coeffs)) deallocate(self % coeffs)

    self % coeffs = coeffs ! implicit allocation
  end subroutine
    
  !!
  !! Calculate release at energy E_in
  !!
  pure function releaseAt(self,E_in) result(release)
    class(polynomialRelease), intent(in) :: self
    real(defReal), intent(in)            :: E_in
    real(defReal)                        :: release
    integer(shortInt)                    :: i

    ! Horner Method Evaluation of a Polynomial
    release = 0.0
    do i=size(self % coeffs),1,-1
      release = release * E_in + self % coeffs(i)
    end do

  end function releaseAt

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(polynomialRelease), intent(inout) :: self

    if(allocated(self % coeffs)) deallocate(self % coeffs)

  end subroutine kill

  !!
  !! Constructor
  !!
  function new_polynomialRelease(coeffs) result(new)
    real(defReal),dimension(:),intent(in)   :: coeffs
    type(polynomialRelease)                 :: new

    call new % init(coeffs)

  end function new_polynomialRelease

  !!
  !! Constructor from ACE
  !! Head of aceCard needs to be set to the beginning of the data (KNU+1 in Table F-5 of MCNP Manual)
  !! NOTE : Defining another init for ACE would help to avoid unnecesary reallocation of memory
  !!
  function new_polynomialRelease_fromACE(ACE) result(new)
    type(aceCard), intent(inout)           :: ACE
    type(polynomialRelease)                :: new
    real(defReal),dimension(:),allocatable :: coeffs
    integer(shortInt)                      :: N

    ! Read number of coefficients
    N = ACE % readInt()

    ! Read coefficients
    coeffs = ACE % readRealArray(N)

    ! Initialise
    call new % init(coeffs)

  end function new_polynomialRelease_fromACE


end module polynomialRelease_class
