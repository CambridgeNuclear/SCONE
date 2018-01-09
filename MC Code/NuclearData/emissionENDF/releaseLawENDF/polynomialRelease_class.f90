module polynomialrelease_class

  use numPrecision
  use releaseLawENDF_class, only : releaseLawENDF

  implicit none
  private

  interface polynomialRelease
    module procedure new_polynomialRelease
  end interface

  type, public,extends(releaseLawENDF) :: polynomialRelease
      private
      real(defReal),dimension(:),allocatable :: coefficients  !! Polynomial coefficients [a_0,a_1,...]
    contains
      procedure :: init
      procedure :: releaseAt
  end type polynomialRelease

contains

  subroutine init(self,coefficients)
    class(polynomialRelease), intent(inout)  :: self
    real(defReal),dimension(:),intent(in)    :: coefficients

    if (allocated(self % coefficients)) deallocate(self % coefficients)

    self % coefficients = coefficients ! implicit allocation
  end subroutine
    
  pure function releaseAt(self,energy) result(release)
    class(polynomialRelease), intent(in)  :: self
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


  function new_polynomialRelease(coefficients) result (newPolynomialRelease)
    real(defReal),dimension(:),intent(in)   :: coefficients
    type(polynomialRelease), pointer         :: newPolynomialRelease

    allocate(newPolynomialRelease)
    call newPolynomialRelease % init(coefficients)

  end function new_polynomialRelease


end module polynomialRelease_class
