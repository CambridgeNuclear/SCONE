module xsCDF_class

  use numPrecision

  implicit none
  private

  type,abstract, public :: xsCDF
    private
  contains
    procedure(invert),deferred  :: invert
  end type xsCDF

  abstract interface

    function invert(self,r) result(mask)
      import :: xsCDF,   &
                defReal, &
                shortInt
      class(xsCDF), intent(in)  :: self
      real(defReal),intent(in)  :: r
      integer(shortInt)         :: mask
    end function

  end interface


contains


    
end module xsCDF_class
