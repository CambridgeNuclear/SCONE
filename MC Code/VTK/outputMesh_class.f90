!!
!! Module used to output mesh data in VTK format
!!
module outputMesh_class

  use numPrecision
  use universalVariables
  use genericProcedures
  use outputVTK_class

  implicit none
  private

  !!
  !! Object responsible for outputting VTK files
  !!
  type, extends(outputVTK), public                 :: outputMesh
    real(defReal), dimension(3)                    :: corner  ! corner of the mesh
    real(defReal), dimension(3)                    :: width   ! mesh cell width in each direction
    integer(shortInt), dimension(3)                :: nVox    ! number of voxels
    integer(shortInt)                              :: nCells  ! total number of mesh cells
    integer(shortInt)                              :: nOutput ! number of separate mesh data outputs
    real(defReal), dimension(:,:,:,:), allocatable :: values  ! mesh values indexed by output number and mesh indices
    character(nameLen), dimension(:), allocatable  :: dataName! name of the dataset
  contains
    procedure :: init
    procedure :: addData => addData
    procedure :: output => output
  end type

contains

  !!
  !! Initialise mesh output by providing mesh structure
  !!
  subroutine init(self, corner, width, nVox)
    class(outputMesh), intent(inout)            :: self
    real(defReal), dimension(3), intent(in)     :: corner
    real(defReal), dimension(3), intent(in)     :: width
    integer(shortInt), dimension(3), intent(in) :: nVox

    self % legacy = .TRUE.
    self % corner = corner
    self % width = width
    self % nVox = nVox
    self % nCells = nVox(1) * nVox(2) * nVox(3)
    self % nOutput = 0

  end subroutine init

  !!
  !! Add a new mesh data set to the object
  !!
  subroutine addData(self, newValues, newName)
    class(outputMesh), intent(inout)               :: self
    real(defReal), dimension(:,:,:), intent(in)    :: newValues
    character(*), intent(in)                       :: newName
    real(defReal), dimension(:,:,:,:), allocatable :: tempArray
    character(nameLen), dimension(:), allocatable  :: tempChar


    ! Perform error checks on the size of the array provided
    if ((size(newValues,1) /= self % nVox(1)) .OR. (size(newValues,2) /= self % nVox(2)) .OR. &
    (size(newValues,3) /= self % nVox(3))) call fatalError&
    ('addData, outputMesh','Input array size does not agree with anticipated size')

    self % nOutput = self % nOutput + 1

    ! If this is the first data set, simply copy the values in and be done
    if (self % nOutput == 1) then
      allocate(self % values(1, self % nVox(1), self % nVox(2), self % nVox(3)))
      self % values(1,:,:,:) = newValues
      allocate(self % dataName(1))
      self % dataName(1) = newName
    ! Otherwise, move previously allocated values into a temporary array
    else
      allocate(tempArray(self % nOutput, self % nVox(1), self % nVox(2), self % nVox(3)))
      tempArray(1:self % nOutput - 1, :, :, :) = self % values
      tempArray(self % nOutput, :, :, :) = newValues
      call move_alloc(tempArray, self % values)
      allocate(tempChar(self % nOutput))
      tempChar(1:self % nOutput - 1) = self % dataName
      tempChar(self % nOutput) = newName
      call move_alloc(tempChar, self % dataName)
    end if

  end subroutine addData

  !!
  !! Output all mesh data
  !!
  subroutine output(self, name)
    class(outputMesh), intent(in)                   :: self
    character(*), intent(in)                        :: name
    integer(shortInt)                               :: i, j, k, l, stat
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

    ! Output header - string terminated by character (256 characters maximum)
    write(10,'(A)') 'Voxel output file generated from the SCONE Monte Carlo code'

    ! Output file format - either ASCII or BINARY
    write(10,'(A)') 'ASCII'

    ! Output dataset structure - describe geometry and topology of the dataset
    ! Here this follows the Structured Points standard
    ! Specify number of elements in each direction, specify an origin,
    ! and specify spacing in each direction
    write(10,'(A)') 'DATASET STRUCTURED_POINTS'
    write(10,'(A,A,A,A,A,A)') 'DIMENSIONS ',trim(numToChar(self % nVox(1))),' '&
                       ,trim(numToChar(self % nVox(2))),' ',trim(numToChar(self % nVox(3)))

    write(10,'(A,A,A,A,A,A)') 'ORIGIN ',trim(numToChar(self % corner(1))),' ',&
                              trim(numToChar(self % corner(2))),' ',trim(numToChar(self % corner(3)))

    write(10,'(A,A,A,A,A,A)') 'SPACING ',trim(numToChar(self % width(1))),' ',&
                              trim(numToChar(self % width(2))),' ',trim(numToChar(self % width(3)))

    write(10,'(A,A)') 'POINT_DATA ',numToChar(self % nCells)

    ! Output dataset attributes - begins with POINT_DATA or CELL_DATA followed by number of cells/points
    ! Other keyword/data combinations then define the dataset attribute values (scalar, vectors, tensors...)
    ! This can be modified to output RGB vectors or scalars rather than integers

    ! Loop over each data set

    do l = 1, self % nOutput
      write(10,'(A,A,A)') 'SCALARS ',trim(self % dataName(l)),' float 1'

      ! Possibility: construct and specify lookup table
      write(10,'(A)') 'LOOKUP_TABLE default'

      do k = 1, self % nVox(3)
        do j = 1, self % nVox(2)
          do i = 1, self % nVox(1)
            write(10,'(A)') numToChar(self % values(l,i,j,k))
          end do
        end do
      end do

    end do

    ! Close the file
    close(10)

  end subroutine output

end module outputMesh_class
