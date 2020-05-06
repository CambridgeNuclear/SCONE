!!
!! Module used to output mesh data in VTK format
!!
module outputVTK_class

  use numPrecision
  use universalVariables
  use genericProcedures
  use dictionary_class, only: dictionary

  implicit none
  private

  !!
  !! Object responsible for creating and outputting VTK files
  !!
  type, public                                     :: outputVTK
    logical(defBool)                               :: legacy = .TRUE. ! Is it legacy VTK?
    integer(shortInt), dimension(2)                :: version = [3,0] ! VTK version
    real(defReal), dimension(3)                    :: corner          ! corner of the mesh
    real(defReal), dimension(3)                    :: width           ! mesh cell width in each direction
    integer(shortInt), dimension(3)                :: nVox            ! number of voxels
    integer(shortInt)                              :: nCells          ! total number of mesh cells
    integer(shortInt)                              :: nOutput         ! number of separate mesh data outputs
    real(defReal), dimension(:,:,:,:), allocatable :: values          ! mesh values indexed by output number and mesh indices
    character(nameLen), dimension(:), allocatable  :: dataName        ! name of the dataset
    logical(defBool), dimension(:), allocatable    :: dataReal        ! is data real?(T) or int?(F)
  contains
    procedure :: init
    generic   :: addData => addDataInt,&
                            addDataReal
    procedure :: addDataInt
    procedure :: addDataReal
    procedure :: output
    procedure :: kill
  end type

contains

  !!
  !! Initialise mesh output by providing mesh structure
  !!
  subroutine init(self, vtkDict)
    class(outputVTK), intent(inout)              :: self
    class(dictionary), intent(in)                :: vtkDict
    real(defReal), dimension(:), allocatable     :: corner
    real(defReal), dimension(:), allocatable     :: width
    integer(shortInt), dimension(:), allocatable :: nVox
    character(nameLen) :: here ='init, outputVTK_class.f90'

    self % legacy = .TRUE.
    call vtkDict % get(corner,'corner')
    self % corner = corner
    call vtkDict % get(width,'width')
    self % width = width
    call vtkDict % get(nVox,'vox')
    self % nVox = nVox
    if (size(self % corner) .NE. 3) then
      call fatalError(here,'Voxel plot requires corner to have 3 values')
    endif
    if (size(self % width) .NE. 3) then
      call fatalError(here,'Voxel plot requires width to have 3 values')
    endif
    if (size(self % nVox) .NE. 3) then
      call fatalError(here,'Voxel plot requires vox to have 3 values')
    endif
    self % nCells = self % nVox(1) * self % nVox(2) * self % nVox(3)
    self % nOutput = 0

  end subroutine init

  !!
  !! Add a new mesh data set to the object
  !!
  subroutine addDataReal(self, newValues, newName)
    class(outputVTK), intent(inout)                :: self
    real(defReal), dimension(:,:,:), intent(in)    :: newValues
    character(*), intent(in)                       :: newName
    real(defReal), dimension(:,:,:,:), allocatable :: tempArray
    character(nameLen), dimension(:), allocatable  :: tempChar
    logical(defBool), dimension(:), allocatable    :: tempBool


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
      allocate(self % dataReal(1))
      self % dataReal(1) = .TRUE.
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
      allocate(tempBool(self % nOutput))
      tempBool(1:self % nOutput - 1) = self % dataReal
      tempBool(self % nOutput) = .TRUE.
      call move_alloc(tempBool, self % dataReal)
    end if

  end subroutine addDataReal


  !!
  !! Add a new mesh data set to the object
  !!
  subroutine addDataInt(self, newValues, newName)
    class(outputVTK), intent(inout)                 :: self
    integer(shortInt), dimension(:,:,:), intent(in) :: newValues
    character(*), intent(in)                        :: newName
    real(defReal), dimension(:,:,:,:), allocatable  :: tempArray
    character(nameLen), dimension(:), allocatable   :: tempChar
    logical(defBool), dimension(:), allocatable     :: tempBool


    ! Perform error checks on the size of the array provided
    if ((size(newValues,1) /= self % nVox(1)) .OR. (size(newValues,2) /= self % nVox(2)) .OR. &
    (size(newValues,3) /= self % nVox(3))) call fatalError&
    ('addData, outputMesh','Input array size does not agree with anticipated size')

    self % nOutput = self % nOutput + 1

    ! If this is the first data set, simply copy the values in and be done
    if (self % nOutput == 1) then
      allocate(self % values(1, self % nVox(1), self % nVox(2), self % nVox(3)))
      self % values(1,:,:,:) = real(newValues,defReal)
      allocate(self % dataName(1))
      self % dataName(1) = newName
      allocate(self % dataReal(1))
      self % dataReal(1) = .FALSE.
    ! Otherwise, move previously allocated values into a temporary array
    else
      allocate(tempArray(self % nOutput, self % nVox(1), self % nVox(2), self % nVox(3)))
      tempArray(1:self % nOutput - 1, :, :, :) = self % values
      tempArray(self % nOutput, :, :, :) = real(newValues,defReal)
      call move_alloc(tempArray, self % values)
      allocate(tempChar(self % nOutput))
      tempChar(1:self % nOutput - 1) = self % dataName
      tempChar(self % nOutput) = newName
      call move_alloc(tempChar, self % dataName)
      allocate(tempBool(self % nOutput))
      tempBool(1:self % nOutput - 1) = self % dataReal
      tempBool(self % nOutput) = .FALSE.
      call move_alloc(tempBool, self % dataReal)
    end if

  end subroutine addDataInt

  !!
  !! Output all mesh data
  !!
  subroutine output(self, name)
    class(outputVTK), intent(in) :: self
    character(*), intent(in)     :: name
    integer(shortInt)            :: i, j, k, l, stat
    character(256)               :: filename

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
            if (self % dataReal(l)) then
              write(10,'(A)') numToChar(self % values(l,i,j,k))
            else
              write(10,'(A)') numToChar(int(self % values(l,i,j,k),shortInt))
            endif
          end do
        end do
      end do

    end do

    ! Close the file
    close(10)

  end subroutine output

  subroutine kill(self)
    class(outputVTK), intent(inout) :: self

    deallocate(self % values)
    deallocate(self % dataName)
    self % nCells = 0
    self % nOutput = 0

  end subroutine kill

end module outputVTK_class
