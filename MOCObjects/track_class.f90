module track_class

  use numPrecision
  use universalVariables
  use genericProcedures

  implicit none
  private

  type, public :: track2D
    type(segment), dimension(:), allocatable :: segments
    integer(shortInt) :: idxTrack
    integer(shortInt) :: idxAz
    integer(shortInt) :: nSegs
    real(defReal), dimension(3) :: start
    real(defReal), dimension(3) :: end
  contains
    !! Public Interface
    procedure            :: init
    procedure            :: addSegment
  end type track2D

contains

  subroutine init(self,idxTrack,idxAz)
    class(track2D), intent(inout)    :: self
    integer(shortInt), intent(inout) :: idxTrack
    real(defReal), intent(in)        :: idxAz
    self % idxTrack = idxTrack + 1
    self % idxAz = idxAz
    self % nSegs = 0

  end subroutine init

  subroutine addSegment(self,seg)
    class(track2D), intent(inout) :: self
    class(segment), intent(in)    :: seg

    self % nSegs = self % nSegs + 1
  end subroutine addSegment

end module track_class

