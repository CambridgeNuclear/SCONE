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
  type, abstract, public :: outputVTK
    logical(defBool)                :: legacy = .TRUE. ! Determines if output is in legacy format
    integer(shortInt), dimension(2) :: version = [3,0] ! VTK
    logical(defBool)                :: initialise = .FALSE.
  contains
    procedure(addData), deferred :: addData
    procedure(output), deferred  :: output
  end type

  abstract interface

    subroutine addData(self, newValues, newName)
      use numPrecision
      import :: outputVTK
      implicit none
      class(outputVTK), intent(inout)             :: self
      real(defReal), dimension(:,:,:), intent(in) :: newValues
      character(nameLen), intent(in)              :: newName
    end subroutine addData

    subroutine output(self, name)
      use numPrecision
      import :: outputVTK
      implicit none
      class(outputVTK), intent(in)   :: self
      character(nameLen), intent(in) :: name
    end subroutine output

  end interface

contains


!  !!
!  !! Output polyline data in VTK legacy format
!  !! Used for, e.g., particle tracks
!  !!
!  subroutine outputPolyLine(self, name, polydata, dataName)
!    class(outputVTK), intent(in)              :: self
!    character(*), intent(in)                  :: name
!    real(defReal), dimension(:,:), intent(in) :: polydata
!    character(*), intent(in), optional        :: dataName
!    integer(shortInt)                         :: nPoints, i, stat
!    character(256)                            :: filename
!
!    nPoints = size(polyData,1)
!    if (size(polyData,2) /= 3) call fatalError&
!    ('outputPolyLine','Line data is not correctly dimensioned')
!
!    ! Append .vtk to filename
!    filename = trim(name)//'.vtk'
!
!    ! Check if file already exists - if so, delete it
!    open(unit = 10, iostat = stat, file = filename, status = 'old')
!    if (stat == 0) close(10, status = 'delete')
!
!    ! Open file to begin output
!    open(unit = 10, file = filename, status = 'new')
!
!    ! Output file version and identifier  ' # vtk DataFile Version x.x
!    write(10,'(A,I0,A,I0)') '# vtk DataFile Version ',self %version(1),'.', self % version(2)
!
!    ! Output header - string terminated by character (256 characters maximum)
!    write(10,'(A)') 'PolyData output file generated from the SCONE Monte Carlo code'
!
!    ! Output file format - either ASCII or BINARY
!    write(10,'(A)') 'ASCII'
!
!    ! Output dataset structure - describe geometry and topology of the dataset
!    ! Here this follows the PolyData standard - specify points and lines
!    ! and specify spacing in each direction
!    write(10,'(A)') 'DATASET POLYDATA'
!    write(10,'(A,A,A)') 'POINTS ',numToChar(nPoints),' float'
!    ! Provide each of the points composing the lines
!    do i = 1,nPoints
!      write(10,'(A,A,A,A,A)') numToChar(polydata(i,1)),' ',numToChar(polydata(i,2)),' ',&
!                              numToChar(polydata(i,3))
!    end do
!    write(10,'(A,A,A)') 'LINES ',numToChar(nPoints - 1),' 1'
!    ! Construct lines from the points
!    do i = 1,nPoints - 1
!      !write(10,'(A,A,A,A)') '2 ',
!    end do
!
!  end subroutine outputPolyLine


end module outputVTK_class

