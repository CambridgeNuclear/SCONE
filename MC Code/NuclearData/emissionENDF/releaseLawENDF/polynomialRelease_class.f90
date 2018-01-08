module polynomialrelease_class

  use numPrecision
  use releaseLawENDF_class, only : releaseLawENDF

  implicit none
  private

  interface polynomialrelease
    module procedure new_polynomialrelease
  end interface

  type, public,extends(releaseLawENDF) :: polynomialrelease
      private
      real(defReal),dimension(:),allocatable :: coefficients  !! Polynomial coefficients [a_0,a_1,...]
    contains
      procedure :: init
      procedure :: releaseAt
  end type polynomialrelease

contains

  subroutine init(self,coefficients)
    class(polynomialrelease), intent(inout)  :: self
    real(defReal),dimension(:),intent(in)   :: coefficients

    if (allocated(self % coefficients)) deallocate(self % coefficients)

    self % coefficients = coefficients ! implicit allocation
  end subroutine
    
  pure function releaseAt(self,energy) result(release)
    class(polynomialrelease), intent(in)  :: self
    real(defReal), intent(in)            :: energy
    real(defReal)                        :: release
    real(defReal)                        :: x
    integer(shortInt)                    :: i


    ! Horner Method Evaluation of a Polynomial
    release = 0.0
    do i=size(self % coefficients),1,-1
      release = release * energy + self % coefficients(i)
    end do

  end function releaseAt


  function new_polynomialrelease(coefficients) result (newPolynomialrelease)
    real(defReal),dimension(:),intent(in)   :: coefficients
    type(polynomialrelease), pointer         :: newPolynomialrelease

    allocate(newPolynomialrelease)
    call newPolynomialrelease % init(coefficients)

  end function new_polynomialrelease


end module polynomialrelease_class
