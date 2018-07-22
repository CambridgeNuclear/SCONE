module fillArray_class

  use numPrecision

  implicit none
  private


  type, public :: fillArray
    private
    integer(shortInt),dimension(:),allocatable :: fills
  contains
    ! Build procedure
    procedure :: init
    procedure :: setFill
    procedure :: getFill
  end type fillArray

contains

  subroutine init(self)
    class(fillArray), intent(inout) :: self
  end subroutine init

  subroutine setFill(self)
    class(fillArray), intent(inout) :: self
  end subroutine

  function getFill(self,uniqueId) result(fill)
    class(fillArray), intent(inout) :: self
    integer(shortInt)               :: uniqueID
  end function getFill


end module fillArray_class
