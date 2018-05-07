!!
!! Object for storing and handling particle tracking information
!! Used for visualisation and bug-checking purposes
!!
module tracking_class
  use genericProcedures
  use numPrecision
  use universalVariables

  use geometry_class

  implicit none

  type, public :: tracking
    real(defReal), private, dimension(3,:) :: points    ! Must modify, should have a linked list format
    integer(shortInt), private, dimension(:) :: discontinuities
  contains
    procedure :: addPoint
    procedure :: addDiscontinuity ! For periodic translations
    procedure :: outputTracks
  end type tracking

contains

  !!
  !! Adds another point to the point array
  !!
  subroutine addPoint(self, r)
    class(tracking), intent(inout) :: self
    real(defReal), dimension(3), intent(in) :: r

  end subroutine addPoint

  !!
  !! Adds an index to the discontinuity array when a translation occurs
  !!
  subroutine addDiscontinuity(self)
    class(tracking), intent(inout) :: self

  end subroutine addDiscontinuity

  !!
  !! Outputs tracking information
  !!
  subroutine outputTracks
    class(tracking), intent(inout) :: self

  end subroutine outputTracks

end module tracking_class
