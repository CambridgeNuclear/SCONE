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
  !! Output voxel information in VTK legacy format
  !!
  subroutine outputVoxels(self, name, arraySize, values, corner, width, dataName)
    class(outputVTK), intent(in)                    :: self
    character(*), intent(in)                        :: name
    integer(shortInt), dimension(3), intent(in)     :: arraySize
    integer(shortInt), dimension(:,:,:), intent(in) :: values
    real(defReal), dimension(3), intent(in)         :: corner
    real(defReal), dimension(3), intent(in)         :: width
    character(*), intent(in), optional              :: dataName ! Optional identifier for the data - mustn't have whitespace!!!!
    integer(shortInt)                               :: i, j, k, stat
    character(256)                                  :: filename

    ! Append .vtk to filename
    filename = trim(name)//'.vtk'

    ! Check if file already exists - if so, delete it
    open(unit = 10, iostat = stat, file = filename, status = 'old')
    if (stat == 0) close(10, status = 'delete')

    ! Open file to begin output
    open(unit = 10, file = filename, status = 'new')

    ! Output file version and identifier  ' # vtk DataFile Version x.x
    write(10,'(A,I0,A,I0)') '# vtk DataFile Version ',self %version(1),'.', self % version(2)

    ! Output header - string terminated by character \n (256 characters maximum)
    write(10,'(A)') 'Voxel output file generated from the SCONE Monte Carlo code \n'

    ! Output file format - either ASCII or BINARY
    write(10,'(A)') 'ASCII'

    ! Output dataset structure - describe geometry and topology of the dataset
    ! Begins with a line containing DATASET followed by a keyword for the dataset type
    ! Other keyword/data combindations define the actual data
    !
    ! Here this follows the Structured Points standard
    ! Specify number of elements in each direction, specify an origin,
    ! and specify spacing in each direction
    write(10,'(A)') 'DATASET STRUCTURED_POINTS'
    write(10,'(A,A,A,A,A,A)') 'DIMENSIONS ',trim(numToChar(arraySize(1))),' '&
                       ,trim(numToChar(arraySize(2))),' ',trim(numToChar(arraySize(3)))

    write(10,'(A,A,A,A,A,A)') 'ORIGIN ',trim(numToChar(corner(1))),' ',trim(numToChar(corner(2)))&
                   ,' ',trim(numToChar(corner(3)))

    write(10,'(A,A,A,A,A,A)') 'SPACING ',trim(numToChar(width(1))),' ',trim(numToChar(width(2))),&
                              ' ',trim(numToChar(width(3)))

    ! Output dataset attributes - begins with POINT_DATA or CELL_DATA followed by number of cells/points
    ! Other keyword/data combinations then define the dataset attribute values (scalar, vectors, tensors...)
    ! This can be modified to output RGB vectors or scalars rather than integers

    ! If dataName is present, apply the name
    ! Also, assume output is an integer array (for now!)
    if(present(dataName)) then
      write(10,'(A,A,A)') 'SCALARS ',trim(dataName),' int 1'
    else
      write(10,'(A)') 'SCALARS int 1'
    end if

    ! Possibility: construct and specify lookup table
    write(10,'(A)') 'LOOKUP_TABLE default'

    do k = 1, arraySize(3)
      do j = 1, arraySize(2)
        do i = 1, arraySize(1)
          write(10,'(I0)') values(i,j,k)
        end do
      end do
    end do

    ! Close the file
    close(10)

  end subroutine outputVoxels


end module outputVTK_class

