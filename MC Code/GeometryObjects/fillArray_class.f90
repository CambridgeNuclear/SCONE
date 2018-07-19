module fillArray_class

  use numPrecision

  implicit none
  private

  integer(shortInt), parameter :: UNSET = 0

  type,private :: fillContent
    integer(shortInt) :: fill = 0
    integer(shortInt) :: ptr  = UNSET

  end type fillContent


  type, public :: fillArray
    private
    type(fillContent), dimension(:), allocatable :: fills
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
