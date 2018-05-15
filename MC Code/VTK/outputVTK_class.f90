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
    integer(shortInt)                               :: i, j, k
    character(256)                                  :: filename, &
                                                       versionAndIdentifier, &
                                                       outputHeader, &
                                                       fileFormat, &
                                                       dataStructure, &
                                                       outputOrigin, &
                                                       outputSpacing, &
                                                       outputDimensions, &
                                                       datasetAttributes, &
                                                       lookupTable

    ! Append .vtk to filename
    filename = trim(name)//'.vtk'

    ! Open file to begin output
    open(unit = 10, file = filename, status = 'new')

    ! Output file version and identifier  ' # vtk DataFile Version x.x
    versionAndIdentifier = '# vtk DataFile Version '//numToChar(self %version(1))//'.'//&
                           numToChar(self % version(2))
    write(10,*) trim(versionAndIdentifier)

    ! Output header - string terminated by character \n (256 characters maximum)
    outputHeader = 'Voxel output file generated from the SCONE Monte Carlo code \n'
    write(10,*) trim(outputHeader)

    ! Output file format - either ASCII or BINARY
    fileFormat = 'ASCII'
    write(10,*) trim(fileFormat)

    ! Output dataset structure - describe geometry and topology of the dataset
    ! Begins with a line containing DATASET followed by a keyword for the dataset type
    ! Other keyword/data combindations define the actual data
    !
    ! Here this follows the Structured Points standard
    ! Specify number of elements in each direction, specify an origin,
    ! and specify spacing in each direction
    dataStructure = 'DATASET STRUCTURED_POINTS'
    outputDimensions = 'DIMENSIONS '//numToChar(arraySize(1))//' '//numToChar(arraySize(2))//&
                       ' '//numToChar(arraySize(3))
    write(10,*) trim(outputDimensions)

    outputOrigin = 'ORIGIN '//numToChar(corner(1))//' '//numToChar(corner(2))&
                   //' '//numToChar(corner(3))
    write(10,*) trim(outputOrigin)
    outputSpacing = 'SPACING '//numToChar(width(1))//' '//numToChar(width(2))//' '//numToChar(width(3))
    write(10,*) trim(outputSpacing)

    ! Output dataset attributes - begins with POINT_DATA or CELL_DATA followed by number of cells/points
    ! Other keyword/data combinations then define the dataset attribute values (scalar, vectors, tensors...)
    ! This can be modified to output RGB vectors or scalars rather than integers

    ! If dataName is present, apply the name
    ! Also, assume output is an integer array (for now!)
    if(present(dataName)) then
      datasetAttributes = 'SCALARS '//trim(dataName)//' int 1'
    else
      datasetAttributes = 'SCALARS int 1'
    end if
    write(10,*) trim(datasetAttributes)

    ! Possibility: construct and specify lookup table
    lookupTable = 'LOOKUP_TABLE default'
    write(10,*) trim(lookupTable)

    do k = 1, arraySize(3)
      do j = 1, arraySize(2)
        do i = 1, arraySize(1)
          write(10,*) numToChar(values(i,j,k))
        end do
      end do
    end do

    ! Close the file
    close(10)

  end subroutine outputVoxels


end module outputVTK_class

