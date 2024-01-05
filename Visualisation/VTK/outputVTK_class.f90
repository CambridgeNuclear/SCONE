!!
!! Module used to output mesh data in VTK format
!!
module outputVTK_class

  use numPrecision
  use universalVariables
  use genericProcedures, only: fatalError
  use dictionary_class, only: dictionary

  implicit none
  private

  !!
  !! Object responsible for creating and outputting VTK files
  !!
  !! Receives data sets and geometric information, generating a legacy VTK file
  !!
  !! Private members:
  !!   legacy   -> is it legacy VTK? May be used in future for different output format
  !!   version  -> describes VTK version
  !!   corner   -> corner of the geometry
  !!   width    -> width of mesh in each direction
  !!   nVox     -> number of voxels in each direction
  !!   nCells   -> total number of voxels
  !!   nOutput  -> number of separate outputs
  !!   values   -> values contained in each voxel for each data set
  !!   dataName -> name of each data set
  !!   dataReal -> does a data set contain reals?(T) or ints?(F)
  !!
  !! Interface:
  !!   init        -> initialises outputVTK
  !!   addData     -> adds data to outputVTK, either real or int data
  !!   addDataInt  -> adds integer data
  !!   addDataReal -> adds real data
  !!   output      -> outputs all data to a .vtk file
  !!   kill        -> cleans up outputVTK
  !!
  !! Sample dictionary input:
  !!   myVTKOutput{
  !!     type vtk;
  !!     corner (-1.2 -3.4 5.7);
  !!     width (0.8 27 100.1);
  !!     vox (1000 2000 5000)
  !!     #what material;#
  !!   }
  !!
  type, public                                              :: outputVTK
    logical(defBool), private                               :: legacy = .TRUE. 
    integer(shortInt), dimension(2), private                :: version = [3,0] 
    real(defReal), dimension(3), private                    :: corner 
    real(defReal), dimension(3), private                    :: width  
    integer(shortInt), dimension(3), private                :: nVox   
    integer(shortInt), private                              :: nCells 
    integer(shortInt), private                              :: nOutput
    real(defReal), dimension(:,:,:,:), allocatable, private :: values  
    character(nameLen), dimension(:), allocatable, private  :: dataName
    logical(defBool), dimension(:), allocatable, private    :: dataReal
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
  !! Provides outputVTK with necessary info regarding voxel geometry
  !! prior to providing any data sets
  !!
  !! Args:
  !!   vtkDict [in] -> dictionary containing VTK file details
  !!
  !! Result:
  !!   Initialised outputVTK, ready to receive data
  !!
  !! Errors:
  !!   Will produce errors if geometry descriptors are not three dimensional
  !!   or if the number of voxels are negative
  !!
  subroutine init(self, vtkDict)
    class(outputVTK), intent(inout)              :: self
    class(dictionary), intent(in)                :: vtkDict
    real(defReal), dimension(:), allocatable     :: corner
    real(defReal), dimension(:), allocatable     :: width
    integer(shortInt), dimension(:), allocatable :: nVox
    character(100) :: here ='init (outputVTK_class.f90)'

    self % legacy = .true.
    call vtkDict % get(corner,'corner')
    self % corner = corner
    call vtkDict % get(width,'width')
    self % width = width
    call vtkDict % get(nVox,'vox')
    self % nVox = nVox
    if (size(self % corner) .NE. 3) then
      call fatalError(here,'Voxel plot requires corner to have 3 values')
    end if
    if (size(self % width) .NE. 3) then
      call fatalError(here,'Voxel plot requires width to have 3 values')
    end if
    if (size(self % nVox) .NE. 3) then
      call fatalError(here,'Voxel plot requires vox to have 3 values')
    end if
    if (any(self % nVox < 0)) then
      call fatalError(here,'Number of voxels must be positive')
    end if
    self % nCells = self % nVox(1) * self % nVox(2) * self % nVox(3)
    self % nOutput = 0

  end subroutine init

  !!
  !! Add a new mesh data set to the object
  !!
  !! Adds a new data set of real data to the file to be written
  !!
  !! Args:
  !!   newValues [in] -> real values corresponding to different voxels
  !!   newName [in]   -> name of the data set
  !!
  !! Result:
  !!   New real data set added to VTK file to be written
  !!
  !! Errors:
  !!   Will return an error if the newValues do not have the expected dimensions
  !!
  subroutine addDataReal(self, newValues, newName)
    class(outputVTK), intent(inout)                :: self
    real(defReal), dimension(:,:,:), intent(in)    :: newValues
    character(*), intent(in)                       :: newName
    real(defReal), dimension(:,:,:,:), allocatable :: tempArray
    character(nameLen), dimension(:), allocatable  :: tempChar
    logical(defBool), dimension(:), allocatable    :: tempBool
    character(100) :: here ='addDataReal (outputVTK_class.f90)'

    ! Perform error checks on the size of the array provided
    if( any(shape(newValues) /= self % nVox)) then
     call fatalError&
    (here,'Input array size does not agree with anticipated size')
    end if
    self % nOutput = self % nOutput + 1

    ! If this is the first data set, simply copy the values in and be done
    if (self % nOutput == 1) then
      allocate(self % values(1, self % nVox(1), self % nVox(2), self % nVox(3)))
      self % values(1,:,:,:) = newValues
      allocate(self % dataName(1))
      self % dataName(1) = newName
      allocate(self % dataReal(1))
      self % dataReal(1) = .true.

    else ! move previously allocated values into a temporary array
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
      tempBool(self % nOutput) = .true.
      call move_alloc(tempBool, self % dataReal)
    end if

  end subroutine addDataReal

  !!
  !! Add a new mesh data set to the object
  !!
  !! Adds a new data set of int data to the file to be written
  !!
  !! Args:
  !!   newValues [in] -> int values corresponding to different voxels
  !!   newName [in]   -> name of the data set
  !!
  !! Result:
  !!   New int data set added to VTK file to be written
  !!
  !! Errors:
  !!   Will return an error if the newValues do not have the expected dimensions
  !!
  subroutine addDataInt(self, newValues, newName)
    class(outputVTK), intent(inout)                 :: self
    integer(shortInt), dimension(:,:,:), intent(in) :: newValues
    character(*), intent(in)                        :: newName
    real(defReal), dimension(:,:,:,:), allocatable  :: tempArray
    character(nameLen), dimension(:), allocatable   :: tempChar
    logical(defBool), dimension(:), allocatable     :: tempBool
    character(pathLen) :: here ='addDataInt (outputVTK_class.f90)'

    ! Perform error checks on the size of the array provided
    if( any(shape(newValues) /= self % nVox)) then
     call fatalError&
    (here,'Input array size does not agree with anticipated size')
    end if
    self % nOutput = self % nOutput + 1

    ! If this is the first data set, simply copy the values in and be done
    if (self % nOutput == 1) then
      allocate(self % values(1, self % nVox(1), self % nVox(2), self % nVox(3)))
      self % values(1,:,:,:) = real(newValues, defReal)
      allocate(self % dataName(1))
      self % dataName(1) = newName
      allocate(self % dataReal(1))
      self % dataReal(1) = .false.

    else ! move previously allocated values into a temporary array
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
      tempBool(self % nOutput) = .false.
      call move_alloc(tempBool, self % dataReal)

    end if

  end subroutine addDataInt

  !!
  !! Output all mesh data
  !!
  !! Writes the .vtk file containing all provided data
  !!
  !! Args:
  !!   name [in] -> name of the output file
  !!
  !! Result:
  !!   A .vtk file containing all data contained by the outputVTK object
  !!
  subroutine output(self, name)
    class(outputVTK), intent(in) :: self
    character(*), intent(in)     :: name
    integer(shortInt)            :: l, stat, file
    character(256)               :: filename

    ! Append .vtk to filename
    filename = trim(name)//'.vtk'

    ! Check if file already exists - if so, delete it
    open(newunit = file, iostat = stat, file = filename)
    if (stat == 0) close(file, status = 'delete')

    ! Open file to begin output
    open(newunit = file, file = filename)

    ! Output file version and identifier  ' # vtk DataFile Version x.x
    write(file,'(A,I0,A,I0)') '# vtk DataFile Version ',self %version(1),'.', self % version(2)

    ! Output header - string terminated by character (256 characters maximum)
    write(file,'(A)') 'Voxel output file generated from the SCONE Monte Carlo code'

    ! Output file format - either ASCII or BINARY
    write(file,'(A)') 'ASCII'

    ! Output dataset structure - describe geometry and topology of the dataset
    ! Here this follows the Structured Points standard
    ! Specify number of elements in each direction, specify an origin,
    ! and specify spacing in each direction
    write(file,'(A)') 'DATASET STRUCTURED_POINTS'
    write(file,'(A,I0,A,I0,A,I0)') 'DIMENSIONS ',self % nVox(1),' ',self % nVox(2),' ',self % nVox(3)
    write(file,'(A,F0.3,A,F0.3,A,F0.3)') 'ORIGIN ',self % corner(1),' ',self % corner(2),' ',self % corner(3)
    write(file,'(A,F0.3,A,F0.3,A,F0.3)') 'SPACING ',self % width(1),' ',self % width(2),' ',self % width(3)
    write(file,'(A,I0)') 'POINT_DATA ',self % nCells

    ! Output dataset attributes - begins with POINT_DATA or CELL_DATA followed by number of cells/points
    ! Other keyword/data combinations then define the dataset attribute values (scalar, vectors, tensors...)
    ! This can be modified to output RGB vectors or scalars rather than integers

    ! Loop over each data set

    do l = 1, self % nOutput
      write(file,'(A,A,A)') 'SCALARS ',trim(self % dataName(l)),' float 1'

      ! Possibility: construct and specify lookup table
      write(file,'(A)') 'LOOKUP_TABLE default'

      if (self % dataReal(l)) then
        write(file,'(F0.3)') self % values(l,:,:,:)
      else
        write(file,'(I0)') int(self % values(l,:,:,:),shortInt)
      endif
    end do

    ! Close the file
    close(file)

  end subroutine output

  !!
  !! Kill outputVTK
  !!
  !! Cleans up contents of outputVTK
  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(outputVTK), intent(inout) :: self

    if(allocated(self % values))   deallocate(self % values)
    if(allocated(self % dataName)) deallocate(self % dataName)
    self % nCells = 0
    self % nOutput = 0

  end subroutine kill

end module outputVTK_class
