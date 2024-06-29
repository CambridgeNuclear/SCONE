module visualiser_class

  use numPrecision
  use universalVariables
  use genericProcedures,  only : fatalError, numToChar
  use hashFunctions_func, only : knuthHash, FNV_1
  use imgBmp_func,        only : imgBmp_toFile
  use commandLineUI,      only : getInputFile
  use dictionary_class,   only : dictionary
  use geometry_inter,     only : geometry
  use materialMenu_mod,   only : mm_colourMap => colourMap
  use outputVTK_class

  implicit none
  private

  !!
  !! Object responsible for controlling visualisation
  !!
  !! Object that creates images relating to SCONE geometries
  !! Should be extensible for adding different visualisation methods
  !! Receives and generates data for visualisation
  !! Requires a dictionary input which specifies the procedures to call
  !! Presently supports: VTK voxel mesh creation
  !!
  !! Private members:
  !!   name    -> name to be used for generating output files (corresponding to input)
  !!   geom    -> pointer to geometry
  !!   vizDict -> dictionary containing visualisations to be generated
  !!
  !! Interface:
  !!   init        -> initialises visualiser
  !!   makeViz     -> constructs requested visualisations
  !!   kill        -> cleans up visualiser
  !!   initVTK     -> initialise a VTK file for later data addition
  !!   addVTKData  -> Add a VTK field to a file
  !!   finaliseVTK -> Conclude writing to a VTK and output
  !!
  !! Sample dictionary input:
  !!   viz{
  !!     vizDict1{ <outputVTK> }
  !!     #vizDict2{ <outputVTK> }#
  !!   }
  !!
  !! NOTE: For details regarding contents of the vizDict dictionaries see the documentation
  !!   of 'makeVTK' and 'makeBmpImg' functions
  !!
  type, public :: visualiser
    character(nameLen), private       :: name
    class(geometry), pointer, private :: geom => null()
    type(dictionary), private         :: vizDict
    type(outputVTK), private          :: vtk
  contains
    procedure :: init
    procedure :: makeViz
    procedure :: kill
    procedure, private :: makeVTK
    procedure, private :: makeBmpImg
    procedure :: slicePlot
    procedure :: voxelPlot
    procedure :: initVTK
    procedure :: addVTKData
    procedure :: finaliseVTK
  end type

contains

  !!
  !! Initialises visualiser
  !!
  !! Provides visualiser with filename for output,
  !! geometry information, and the dictionary describing
  !! what is to be plotted
  !!
  !! Args:
  !!   geom [inout] -> pointer to the geometry
  !!   vizDict[in]  -> dictionary containing what is to be visualised
  !!
  !! Result:
  !!   Initialised visualiser
  !!
  subroutine init(self, geom, vizDict, str)
    class(visualiser), intent(inout)         :: self
    class(geometry), pointer, intent(inout)  :: geom
    class(dictionary), intent(in)            :: vizDict
    character(nameLen), intent(in), optional :: str
    character(:), allocatable                :: string

    ! Obtain file name
    if (present(str)) then
      string = str
    else
      call getInputFile(string)
    end if
    self % name = string

    ! Point to geometry
    self % geom => geom

    ! Store visualisation dictionary
    self % vizDict = vizDict

  end subroutine init

  !!
  !! Generate all visualisations specified by vizDict
  !!
  !! Proceed through all dictionaries contained within vizDict
  !! and perform all corresponding visualisations
  !!
  !! Result:
  !!   Visualisation outputs corresponding to dictionary contents
  !!
  !! Errors:
  !!   Returns an error if an unrecognised visualisation is requested
  !!
  subroutine makeViz(self)
    class(visualiser), intent(inout)             :: self
    class(dictionary), pointer                   :: tempDict
    character(nameLen),dimension(:), allocatable :: keysArr
    integer(shortInt)                            :: i
    character(nameLen)                           :: type
    character(nameLen) :: here ='makeViz (visualiser_class.f90)'

    ! Loop through each sub-dictionary and generate visualisation
    ! (if the visualisation method is available)
    call self % vizDict % keys(keysArr,'dict')

    do i=1,size(keysArr)
      tempDict => self % vizDict % getDictPtr(keysArr(i))
      call tempDict % get(type,'type')
      select case(type)
        case('vtk')
          call self % makeVTK(tempDict)

        case('bmp')
          call self % makeBmpImg(tempDict)

        case default
          call fatalError(here, 'Unrecognised visualisation - presently only accept vtk')

      end select

    end do

  end subroutine makeViz

  !!
  !! Generate a VTK output
  !!
  !! Creates the VTK file corresponding to the contents of dict
  !!
  !! Args:
  !!   dict [in] -> dictionary containing description of VTK file to be made
  !!
  !! Sample input dictionary:
  !!   VTK {
  !!     type vtk;
  !!     corner (-1.0 -1.0 -1.0);  // lower corner of the plot volume
  !!     width  (2.0 2.0 2.0);     // width in each direction
  !!     vox (300 300 300);        // Resolution in each direction
  !!     #what uniqueId;#   // Plot target. 'material' or 'uniqueId'. Default: 'material'
  !!   }
  !!
  !! TODO: VTK output is placed in a input filename appended by '.vtk' extension.
  !!   This prevents multiple VTK visualisations (due to overriding). Might also become
  !!   weird for input files with extension e.g. 'input.dat'.
  !!   DEMAND USER TO GIVE OUTPUT NAME
  !!
  subroutine makeVTK(self, dict)
    class(visualiser), intent(inout)                :: self
    class(dictionary), intent(in)                   :: dict
    type(outputVTK)                                 :: vtk
    integer(shortInt), dimension(:,:,:), allocatable:: voxelMat
    real(defReal), dimension(:), allocatable        :: corner  ! corner of the mesh
    real(defReal), dimension(:), allocatable        :: center  ! center of the mesh
    real(defReal), dimension(:), allocatable        :: width   ! corner of the mesh
    integer(shortInt), dimension(:), allocatable    :: nVox    ! number of mesh voxels
    character(nameLen)                              :: what
    character(nameLen) :: here ='makeVTK (visualiser_class.f90)'

    call vtk % init(dict)

    ! Identify whether plotting 'material' or 'cellID'
    call dict % getOrDefault(what, 'what', 'material')

    ! Obtain geometry data
    call dict % get(corner, 'corner')
    call dict % get(width, 'width')
    call dict % get(nVox, 'vox')

    ! Avoid uninitialised warning
    allocate(center(size(corner)))
    center = corner + width/TWO

    if (size(corner) /= 3) then
      call fatalError(here,'Voxel plot requires corner to have 3 values')
    endif
    if (size(width) /= 3) then
      call fatalError(here,'Voxel plot requires width to have 3 values')
    endif
    if (size(nVox) /= 3) then
      call fatalError(here,'Voxel plot requires vox to have 3 values')
    endif
    allocate(voxelMat(nVox(1), nVox(2), nVox(3)))

    ! Have geometry obtain data
    call self % voxelPlot(voxelMat, center, what, width)

    ! In principle, can add multiple data sets to VTK - not done here yet
    ! VTK data set will use 'what' variable as a name
    call vtk % addData(voxelMat, what)
    call vtk % output(self % name)
    call vtk % kill()

  end subroutine makeVTK

  !!
  !! Output an already constructed VTK and clean up
  !!
  subroutine finaliseVTK(self)
    class(visualiser), intent(inout) :: self

    call self % vtk % output(self % name)
    call self % vtk % kill()

  end subroutine finaliseVTK

  !!
  !! Generate a BMP slice image of the geometry
  !!
  !! Args:
  !!   dict [in] -> Dictionary with settings
  !!
  !! Sample dictionary input:
  !!   bmp_img {
  !!     type bmp;
  !!     #what uniqueID;#      // Target of the plot. 'uniqueId' or 'material'. Default: 'material'
  !!     output img;           // Name of output file without extension
  !!     centre (0.0 0.0 0.0); // Coordinates of the centre of the plot
  !!     axis x;               // Must be 'x', 'y' or 'z'
  !!     res (300 300);        // Resolution of the image
  !!     #offset 978; #        // Parameter to 'randomize' the colour map
  !!     #width (1.0 2.0);#    // Width of the plot from the centre
  !!   }
  !!
  !! NOTE: If 'width' is not given, the plot will extend to the bounds of the geometry.
  !!   This may result in the provided centre being moved to the center of the geometry in the
  !!   plot plane. However, the position on the plot axis will be unchanged.
  !!
  subroutine makeBmpImg(self, dict)
    class(visualiser), intent(inout) :: self
    class(dictionary), intent(in)    :: dict
    real(defReal), dimension(3)      :: centre
    real(defReal), dimension(2)      :: width
    character(1)                     :: dir
    character(nameLen)               :: tempChar
    logical(defBool)                 :: useWidth
    character(nameLen)               :: what, outputFile
    real(defReal), dimension(:), allocatable       :: temp
    integer(shortInt), dimension(:), allocatable   :: tempInt
    integer(shortInt), dimension(:,:), allocatable :: img
    integer(shortInt)                              :: offset
    character(10)                                  :: time
    character(8)                                   :: date
    character(100), parameter :: Here = 'makeBmpImg (visualiser_class.f90)'

    ! Get plot parameters

    ! Identify whether plotting 'material' or 'cellID'
    call dict % getOrDefault(what, 'what', 'material')

    ! Get name of the output file
    call dict % get(outputFile, 'output')
    outputFile = trim(outputFile) // '.bmp'

    ! Central point
    call dict % get(temp, 'centre')

    if (size(temp) /= 3) then
      call fatalError(Here, "'centre' must have size 3. Has: "//numToChar(size(temp)))
    end if

    centre = temp

    ! Axis
    call dict % get(tempChar, 'axis')

    if (len_trim(tempChar) /= 1) then
      call fatalError(Here, "'axis' must be x,y or z. Not: "//tempChar)
    end if

    dir = tempChar(1:1)

    ! Resolution
    call dict % get(tempInt, 'res')

    if (size(tempInt) /= 2) then
      call fatalError(Here, "'res' must have size 2. Has: "//numToChar(size(tempInt)))
    else if (any(tempInt <= 0)) then
      call fatalError(Here, "Resolution must be +ve. There is 0 or -ve entry!")
    end if

    allocate(img(tempInt(1), tempInt(2)))

    ! Optional width
    useWidth = dict % isPresent('width')
    if (useWidth) then
      call dict % get(temp, 'width')

      ! Check for errors
      if (size(temp) /= 2) then
        call fatalError(Here, "'width' must have size 2. Has: "//numToChar((size(temp))))
      else if (any(temp <= ZERO)) then
        call fatalError(Here, "'width' must be +ve. It isn't.")
      end if

      width = temp

    end if

    ! Colourmap offset
    ! If not given select randomly
    if (dict % isPresent('offset')) then
      call dict % get(offset, 'offset')

    else
      call date_and_time(date, time)
      call FNV_1(date // time, offset)

    end if

    ! Get plot
    if (useWidth) then
      call self % slicePlot(img, centre, dir, what, width)
    else
      call self % slicePlot(img, centre, dir, what)
    end if

    ! Translate to an image
    select case (what)
      case ('material')
        img = materialColour(img, offset)

      case ('uniqueID')
        img = uniqueIDColour(img)

      case default
        call fatalError(Here, "Invalid request for plot target. Must be 'material' or 'uniqueID'&
                             & is: "//what)
    end select

    ! Print image
    call imgBmp_toFile(img, outputFile)

  end subroutine makeBmpImg

  !!
  !! Terminates visualiser
  !!
  !! Cleans up remnants of visualiser once it is no longer needed
  !!
  !! Result:
  !!   An empty visualiser object
  !!
  subroutine kill(self)
    class(visualiser), intent(inout) :: self

    self % name =''
    self % geom => null()
    call self % vizDict % kill()

  end subroutine kill


  !!
  !! Convert matIdx to a 24bit colour
  !!
  !! Special materials are associated with special colours:
  !!   OUTSIDE_MAT -> white (#ffffff)
  !!   VOID_MAT    -> black (#000000)
  !!   UNDEF_MAT   -> green (#00ff00)
  !!
  !! Args:
  !!   matIdx [in] -> Value of the material index
  !!   offset [in] -> Offset to be used in the hash function
  !!
  !! Result:
  !!   A 24-bit colour specifying the material
  !!
  elemental function materialColour(matIdx, offset) result(colour)
    integer(shortInt), intent(in) :: matIdx
    integer(shortInt)             :: colour
    integer(shortInt), intent(in) :: offset

    ! Since Knuth hash is cheap we can compute it anyway even if colour from
    ! the map will end up being used
    colour = mm_colourMap % getOrDefault(matIdx, knuthHash(matIdx + offset, 24))

  end function materialColour

  !!
  !! Convert uniqueID to 24bit colour
  !!
  !! An elemental wrapper over Knuth Hash
  !!
  !! We use a hash function to scatter colours across all available.
  !! Knuth multiplicative hash is very good at scattering integer
  !! sequences e.g. {1, 2, 3...}. Thus, it is ideal for a colourmap.
  !!
  !! Args:
  !!   uniqueID [in] -> Value of the uniqueID
  !!
  !! Result:
  !!   A 24-bit colour specifying the uniqueID
  !!
  elemental function uniqueIDColour(uniqueID) result(colour)
    integer(shortInt), intent(in) :: uniqueID
    integer(shortInt)             :: colour

    colour = knuthHash(uniqueID, 24)

  end function uniqueIDColour

  !!
  !! Produce a 2D plot of the geometry
  !!
  !! Resolution is determined by a size of provided output matrix
  !! By default plot plane is normal to z-axis, with width determined by bounds of the
  !! geometry.
  !!
  !! Args:
  !!   img [out]   -> Rank 2 matrix. It is effectively a bitmap image
  !!   centre [in] -> Location of the centre of the image
  !!   dir [in]    -> Axis normal to plot plane. In {'x','y','z'}
  !!   what [in]   -> What to plot 'material' or 'uniqueID'
  !!   width [in]  -> Optional. Width of the plot in both directions. Direction lower in
  !!     sequence {x,y,z} is given first.
  !!
  subroutine slicePlot(self, img, centre, dir, what, width)
    class(visualiser), intent(in)                     :: self
    integer(shortInt), dimension(:,:), intent(out)    :: img
    real(defReal), dimension(3), intent(in)           :: centre
    character(1), intent(in)                          :: dir
    character(*), intent(in)                          :: what
    real(defReal), dimension(2), optional, intent(in) :: width
    real(defReal), dimension(3)     :: low, top , step, point, corner
    real(defReal), dimension(6)     :: aabb
    integer(shortInt), dimension(2) :: plane
    integer(shortInt)               :: ax, i, j, matIdx, uniqueID
    logical(defBool)                :: printMat
    character(100), parameter :: Here = 'slicePlot (visualiser_class.f90)'

    ! Select plane of the plot
    select case (dir)
      case ('x')
        ax = X_AXIS
        plane = [Y_AXIS, Z_AXIS]

      case ('y')
        ax = Y_AXIS
        plane = [X_AXIS, Z_AXIS]

      case ('z')
        ax = Z_AXIS
        plane = [X_AXIS, Y_AXIS]

      case default
        call fatalError(Here, 'Unknown normal axis: '//dir)
        ax = X_AXIS
        plane = 0 ! Make compiler happy
    end select

    ! Find lower and upper corner of the plot
    if (present(width)) then
      low(plane) = centre(plane) - width * HALF
      top(plane) = centre(plane) + width * HALF
      low(ax) = centre(ax)
      top(ax) = centre(ax)

    else
      aabb = self % geom % bounds()
      low = aabb(1:3)
      top = aabb(4:6)
      low(ax) = centre(ax)
      top(ax) = centre(ax)

    end if

    ! Calculate step size in all directions
    step(ax) = ZERO
    step(plane) = (top(plane) - low(plane)) / shape(img)

    ! Select what to print
    select case (what)
      case ('material')
        printMat = .true.

      case('uniqueID')
        printMat = .false.

      case default
        call fatalError(Here, 'Target of plot must be material or uniqueID. Not: '//trim(what))
        printMat = .false. ! Make compiler happy

    end select

    ! Print the image
    corner = low - HALF * step
    point(ax) = corner(ax)

    !$omp parallel do firstprivate(point) private(matIdx, uniqueID)
    do j = 1, size(img, 2)
      point(plane(2)) = corner(plane(2)) + step(plane(2)) * j

      do i = 1, size(img, 1)
        point(plane(1)) = corner(plane(1)) + step(plane(1)) * i

        ! Find material and paint image
        call self % geom % whatIsAt(matIdx, uniqueID, point)

        ! Paint the pixel
        if (printMat) then
          img(i, j) = matIdx
        else
          img(i, j) = uniqueID
        end if

      end do
    end do
    !$omp end parallel do

  end subroutine slicePlot

  !!
  !! Produce a 3D Voxel plot of the geometry
  !!
  !! Resolution is determined by a size of provided output matrix
  !! By default, bounds of the plot correspond to the bounds of the geometry
  !!
  !! Args:
  !!   img [out] -> Rank 3 matrix, It is effectively a 3D bitmap
  !!   centre [in] -> Location of the centre of the image
  !!   what [in]   -> What to plot 'material' or 'uniqueID'
  !!   width [in]  -> Optional. Width of the plot in all directions.
  !!
  subroutine voxelPlot(self, img, centre, what, width)
    class(visualiser), intent(in)                     :: self
    integer(shortInt), dimension(:,:,:), intent(out)  :: img
    real(defReal), dimension(3), intent(in)           :: centre
    character(*), intent(in)                          :: what
    real(defReal), dimension(3), optional, intent(in) :: width
    real(defReal), dimension(3)                       :: width_l, centre_l, point
    real(defReal), dimension(6)                       :: aabb
    real(defReal)                                     :: stepZ
    integer(shortInt)                                 :: i

    ! Get local value of width and centre
    if (present(width)) then
      width_l = width
      centre_l = centre

    else
      aabb = self % geom % bounds()
      centre_l = (aabb(1:3) + aabb(4:6)) * HALF
      width_l  = aabb(4:6) - aabb(1:3)

    end if

    ! Calculate step in z direction
    stepZ = width_l(3) / size(img, 3)

    ! Build voxel plot from multiple slice plots
    point = centre_l
    point(3) = centre_l(3) - width_l(3) * HALF - stepZ * HALF
    do i = 1, size(img, 3)
      point(3) = point(3) + stepZ
      call self % slicePlot(img(:,:,i), point, 'z', what, width_l(1:2))

    end do

  end subroutine voxelPlot
  
  !!
  !! Initialise a VTK output with data addition after
  !!
  !! Creates the VTK file corresponding to the contents of vizDict but does not
  !! output it. Allows more complex fields to be plotted, e.g., based on results.
  !! Accepts the first VTK dictionary it finds.
  !!
  !! VTK dictionary is the standard dictionary used above.
  !!
  !! TODO: VTK output is placed in a input filename appended by '.vtk' extension.
  !!   This prevents multiple VTK visualistions (due to overriding). Might also become
  !!   weird for input files with extension e.g. 'input.dat'.
  !!   DEMAND USER TO GIVE OUTPUT NAME
  !!
  subroutine initVTK(self)
    class(visualiser), intent(inout)                :: self
    character(nameLen),dimension(:), allocatable    :: keysArr
    integer(shortInt)                               :: i
    character(nameLen)                              :: vizType
    logical(defBool)                                :: vtkFound
    class(dictionary), pointer                      :: dict
    integer(shortInt), dimension(:,:,:), allocatable:: voxelMat
    real(defReal), dimension(3)                     :: corner  ! corner of the mesh
    real(defReal), dimension(3)                     :: center  ! center of the mesh
    real(defReal), dimension(3)                     :: width   ! corner of the mesh
    real(defReal), dimension(:), allocatable        :: temp    ! temporary vector
    integer(shortInt), dimension(:), allocatable    :: nVox    ! number of mesh voxels
    character(nameLen)                              :: what
    character(nameLen) :: here ='initVTK (visualiser_class.f90)'
    
    ! Loop through each sub-dictionary and generate visualisation
    ! (if the visualisation method is available)
    call self % vizDict % keys(keysArr,'dict')
    vtkFound = .FALSE.
    do i=1,size(keysArr)
      dict => self % vizDict % getDictPtr(keysArr(i))
      call dict % get(vizType,'type')
      if (vizType == 'vtk') then
        vtkFound = .TRUE.
        exit
      end if
    end do

    if (.NOT. vtkFound) call fatalError(Here,'No VTK data provided in dictionary')

    call self % vtk % init(dict)

    ! Identify whether plotting 'material' or 'cellID'
    call dict % getOrDefault(what, 'what', 'material')

    ! Obtain geometry data
    call dict % get(temp, 'corner')
    if (size(temp) /= 3) then
      call fatalError(Here, "'center' must have size 3. Has: "//numToChar(size(temp)))
    end if
    corner = temp

    call dict % get(temp, 'width')
    if (size(temp) /= 3) then
      call fatalError(Here, "'width' must have size 3. Has: "//numToChar(size(temp)))
    end if
    width  = temp

    center = corner + width/TWO
    call dict % get(nVox, 'vox')

    if (size(corner) /= 3) then
      call fatalError(here,'Voxel plot requires corner to have 3 values')
    end if
    if (size(width) /= 3) then
      call fatalError(here,'Voxel plot requires width to have 3 values')
    end if
    if (size(nVox) /= 3) then
      call fatalError(here,'Voxel plot requires vox to have 3 values')
    end if
    allocate(voxelMat(nVox(1), nVox(2), nVox(3)))

    ! Have geometry obtain data
    call self % voxelPlot(voxelMat, center, what, width)

    ! VTK data set will use 'what' variable as a name
    call self % vtk % addData(voxelMat, what)

  end subroutine initVTK
  
  !!
  !! Add additional data to a VTK file based on an array of values.
  !! The array index will correspond to either the material or uniqueID at
  !! a given position in the geometry.
  !!
  !! Assumes the VTK has already been initialised and uses the first VTK
  !! set of values, i.e., index 1, to check which values to add.
  !!
  subroutine addVTKData(self,dataArray, dataName)
    class(visualiser), intent(inout)             :: self
    real(defReal), dimension(:), intent(in)      :: dataArray
    character(nameLen), intent(in)               :: dataName
    integer(shortInt)                            :: i
    integer(shortInt), save                      :: j, k
    real(defReal), dimension(:,:,:), allocatable :: values
    integer(shortInt), dimension(3)              :: nVox
    !$omp threadprivate(j, k)

    nVox = self % vtk % nVox
    allocate(values(nVox(1),nVox(2),nVox(3)))
    !$omp parallel do schedule(static)
    do i = 1, self % vtk % nVox(1)
      do j = 1, self % vtk % nVox(2)
        do k = 1, self % vtk % nVox(3)
          values(i,j,k) = dataArray(int(self % vtk % values(1,i,j,k)))
        end do
      end do
    end do
    !$omp end parallel do

    call self % vtk % addDataReal(values, dataName)

  end subroutine addVTKData


end module visualiser_class
