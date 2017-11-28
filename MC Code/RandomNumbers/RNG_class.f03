module RNG_class

  use numPrecision

  implicit none
  private

  type, public :: RNG
    private

  contains
    procedure :: get
  end type RNG

contains

  function get(self) result(randomNumber)
    class(RNG), intent(inout)          :: self
    real(kind=defReal)                 :: randomNumber

    randomNumber = rand()

  end function
    
end module RNG_class
