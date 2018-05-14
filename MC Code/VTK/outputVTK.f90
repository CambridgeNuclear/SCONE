!!
!! Module used to output data in a VTK format - initially intended for geometry and track
!! visualisation in ParaView using the legacy VTK format
!!
module outputVTK_class

  use numPrecision
  use universalVariables
  use genericProcedures

  implicit none
  private

  !!
  !! Object responsible for outputting VTK files
  !!
  type, public :: outputVTK
    logical(defBool)                :: legacy = .TRUE. ! Determines if output is in legacy format
    integer(shortInt), dimension(2) :: version = [3,0] ! VTK
  contains
    procedure :: init
    procedure :: outputVoxels
  end type

contains

  !!
  !! Initialise object
  !! Will allow swapping to modern VTK
  !!
  subroutine init(self)
    class(outputVTK), intent(inout) :: self
    self % legacy = .TRUE.
  end subroutine init

  !!
  !! Output voxel information
  !!
  subroutine outputVoxels(self, filename, points, values)

  end subroutine outputVoxels

end module outputVTK_class

