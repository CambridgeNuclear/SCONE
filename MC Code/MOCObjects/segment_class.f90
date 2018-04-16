module segment_class

  use numPrecision
  use universalVariables
  use genericProcedures
  use coord_class

  implicit none
  private

  type, public        :: segment
    real(defReal)     :: length
    integer(shortInt) :: idxFSR
  contains
    procedure         :: build
  end type segment

contains

  subroutine init(self,length,idxFSR)
    class(segment), intent(inout) :: self
    real(defReal), intent(in)     :: length
    real(defReal), intent(in)     :: idxFSR
    self % length = length
    self % cellIdx = idxFSR
  end subroutine init

end module segment_class
