module indexStorageMT_class

  use numPrecision
  use genericProcedures, only : fatalError, linFind, searchError

  implicit none
  private

  !!
  !! Type to store indices of MT reactions in nuclide data
  !! For now the implementation is very crude. Will be improved
  !!
  type, public :: indexStorageMT
    private
    integer(shortInt), dimension(:), allocatable :: MTs
  contains
    procedure :: init
    procedure :: getIdx
  end type indexStorageMT

contains

  !!
  !! Initialise by providing array of MT numbers
  !!
  subroutine init(self,MTs)
    class(indexStorageMT), intent(inout)       :: self
    integer(shortInt),dimension(:), intent(in) :: MTs

    self % MTs = MTs

  end subroutine init

  !!
  !! Returns index of an MT number
  !!
  function getIdx(self,MT) result(idx)
    class(indexStorageMT), intent(in) :: self
    integer(shortInt), intent(in)     :: MT
    integer(shortInt)                 :: idx
    character(100),parameter :: Here='getIdx (indexStorageMT_class.f90)'

    idx = linFind(self % MTs,MT)
    call searchError(idx,Here)
    
  end function getIdx
end module indexStorageMT_class
