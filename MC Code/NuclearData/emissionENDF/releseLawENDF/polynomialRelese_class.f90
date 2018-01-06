module polynomialRelese_class

  use numPrecision
  use releseLawENDF_class, only : releseLawENDF

  implicit none
  private

  interface polynomialRelese
    module procedure new_polynomialRelese
  end interface

  type, public,extends(releseLawENDF) :: polynomialRelese
      private
      real(defReal),dimension(:),allocatable :: coefficients  !! Polynomial coefficients [a_0,a_1,...]
    contains
      procedure :: init
      procedure :: releseAt
  end type polynomialRelese

contains

  subroutine init(self,coefficients)
    class(polynomialRelese), intent(inout)  :: self
    real(defReal),dimension(:),intent(in)   :: coefficients

    if (allocated(self % coefficients)) deallocate(self % coefficients)

    self % coefficients = coefficients ! implicit allocation
  end subroutine
    
  pure function releseAt(self,energy) result(relese)
    class(polynomialRelese), intent(in)  :: self
    real(defReal), intent(in)            :: energy
    real(defReal)                        :: relese
    real(defReal)                        :: x
    integer(shortInt)                    :: i


    ! Horner Method Evaluation of a Polynomial
    relese = self % coefficients(size(self % coefficients))
    do i=size(self % coefficients)-1,1,-1
      relese = relese * energy + self % coefficients(i)
    end do

  end function releseAt


  function new_polynomialRelese(coefficients) result (newPolynomialRelese)
    real(defReal),dimension(:),intent(in)   :: coefficients
    type(polynomialRelese), pointer         :: newPolynomialRelese

    allocate(newPolynomialRelese)
    call newPolynomialRelese % init(coefficients)

  end function new_polynomialRelese


end module polynomialRelese_class
