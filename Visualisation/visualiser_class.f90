module visualiser_class

  use numPrecision
  use universalVariables
  use genericProcedures,  only : fatalError, numToChar, crossProduct
  use hashFunctions_func, only : knuthHash, FNV_1
  use imgBmp_func,        only : imgBmp_toFile
  use commandLineUI,      only : getInputFile
  use dictionary_class,   only : dictionary
  use geometry_inter,     only : geometry
  use materialMenu_mod,   only : mm_colourMap => colourMap, mm_nameMap => nameMap
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
  !!   init    -> initialises visualiser
  !!   makeViz -> constructs requested visualisations
  !!   kill    -> cleans up visualiser
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
  contains
    procedure :: init
    procedure :: makeViz
    procedure :: kill
    procedure, private :: makeVTK
    procedure, private :: makeBmpImg
    procedure, private :: makeRayPlot
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
  subroutine init(self, geom, vizDict)
    class(visualiser), intent(inout)        :: self
    class(geometry), pointer, intent(inout) :: geom
    class(dictionary), intent(in)           :: vizDict
    character(:), allocatable               :: string

    ! Obtain file name
    call getInputFile(string)
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
        
        case('ray')
          call self % makeRayPlot(tempDict)

        case default
          call fatalError(here, 'Unrecognised visualisation - presently only accept vtk, bmp, and ray')

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
    center = corner + width/TWO
    call dict % get(nVox, 'vox')

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
    call self % geom % voxelPlot(voxelMat, center, what, width)

    ! In principle, can add multiple data sets to VTK - not done here yet
    ! VTK data set will use 'what' variable as a name
    call vtk % addData(voxelMat, what)
    call vtk % output(self % name)
    call vtk % kill()

  end subroutine makeVTK

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
      call self % geom % slicePlot(img, centre, dir, what, width)
    else
      call self % geom % slicePlot(img, centre, dir, what)
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
  !! Generate a ray traced image of the geometry using Phong's approximation
  !!
  !! Follows the papers from Gavin Ridley and Brian Nease at SNA+MC 2024.
  !!
  !! The user can input the point *towards which* the camera is pointed, the point from which
  !! the target is viewed, the horizontal field of view (defaults to recommended 70 degrees),
  !! the point from which the light originates, and the number of pixels in 2D.
  !! Additionally, the user may specify which materials are transparent. Otherwise the plot may
  !! be quite boring.
  !!
  !! Optionally, the user can specify which direction is 'up' to determine the orientation of
  !! the camera, the relative amount of diffuse light, the field of view, and the colour offset.
  !!
  !! Args:
  !!   dict [in] -> Dictionary with settings
  !!
  !! Sample dictionary input:
  !!   ray_plot {
  !!     type ray;
  !!     output img;             // Name of output file without extension
  !!     centre (0.0 0.0 0.0);   // Coordinates of the centre/target of the plot
  !!     camera (0.0 0.0 0.0);   // Coordinates of the camera
  !!     res (300 300);          // Resolution of the image
  !!     #light  (0.0 0.0 0.0);# // Coordinates of the light source, defaults to camera position
  !!     #up (0.0 0.0 1.0);#     // Which way is 'up'? Determines the rotation of the camera
  !!     #diffuse 0.3;#          // Fraction of light which is diffuse rather than direct
  !!     #fov 70;#               // Field-of-view in the horizontal axis in degrees
  !!     #offset 978; #          // Parameter to 'randomize' the colour map
  !!     #transparent (mat1, mat2, ... );# // Names of transparent materials
  !!   }
  !!
  !!
  subroutine makeRayPlot(self, dict)
    class(visualiser), intent(inout) :: self
    class(dictionary), intent(in)    :: dict
    real(defReal)                    :: fov, diffuse
    real(defReal), dimension(3) :: centre, up, camera, light, cv, ch, d
    real(defReal), dimension(:), allocatable       :: temp
    integer(shortInt), dimension(:), allocatable   :: tempInt
    character(nameLen), dimension(:), allocatable  :: matNames
    integer(shortInt), dimension(:), allocatable   :: mats
    integer(shortInt), dimension(:,:), allocatable :: img, matIDs
    real(defReal), dimension(:,:), allocatable     :: lum
    real(defReal), dimension(3,3)                  :: M
    integer(shortInt)                              :: i, offset
    character(10)                                  :: time
    character(8)                                   :: date
    character(nameLen)                             :: outputFile
    character(100), parameter :: Here = 'makeRayPlot (visualiser_class.f90)'
    
    ! Get name of the output file
    call dict % get(outputFile, 'output')
    outputFile = trim(outputFile) // '.bmp'

    ! Central/target point
    call dict % get(temp, 'centre')

    if (size(temp) /= 3) then
      call fatalError(Here, "'centre' must have size 3. Has: "//numToChar(size(temp)))
    end if

    centre = temp
    deallocate(temp)

    ! Camera origin
    call dict % get(temp, 'camera')

    if (size(temp) /= 3) then
      call fatalError(Here, "'camera' must have size 3. Has: "//numToChar(size(temp)))
    end if
    camera = temp
    deallocate(temp)
    
    ! Get light location or default to camera location
    if (dict % isPresent('light')) then
      call dict % get(temp, 'light')
    
      if (size(temp) /= 3) then
        call fatalError(Here, "'light' must have size 3. Has: "//numToChar(size(temp)))
      end if
      light = temp
      deallocate(temp)

    else
      light = camera
    end if
    
    ! The up direction, which sets the camera rotation
    if (dict % isPresent('up')) then
      call dict % get(temp, 'up')

      if (size(temp) /= 3) then
        call fatalError(Here, "'up' must have size 3. Has: "//numToChar(size(temp)))
      end if
      up = temp
      up = up / norm2(up)
      deallocate(temp)

    else
      up = [0, 0, 1]
    end if

    ! Optional field of view in degrees
    ! Convert to radians
    call dict % getOrDefault(fov, 'fov', 70.0_defReal)
    fov = fov * PI / 180

    ! Optional fraction of diffuse light
    call dict % getOrDefault(diffuse, 'diffuse', 0.3_defReal)
    if ((diffuse > 1) .or. (diffuse < 0)) call fatalError(Here,'The fraction of diffuse light must be between 0 and 1')

    ! Create governing vectors for the plot
    d = (centre - camera)
    d = d/norm2(d)
    
    ! Ensure that up is not colinear with the view direction
    if (all(abs(crossProduct(up, d)) < 1E-6)) call fatalError(Here,"View direction is co-linear with 'up'.")
    
    cv = crossProduct(d, up)
    cv = cv / norm2(cv)
    
    ch = crossProduct(cv, d)
    ch = ch / norm2(ch)

    ! Create coordinate matrix
    M(:,1) = d
    M(:,2) = cv
    M(:,3) = ch

    ! Resolution
    call dict % get(tempInt, 'res')

    if (size(tempInt) /= 2) then
      call fatalError(Here, "'res' must have size 2. Has: "//numToChar(size(tempInt)))
    else if (any(tempInt <= 0)) then
      call fatalError(Here, "Resolution must be +ve. There is 0 or -ve entry!")
    end if

    allocate(matIDs(tempInt(1), tempInt(2)))
    allocate(lum(tempInt(1), tempInt(2)))

    ! Transparent materials
    ! Defaults to an array containing only VOID and OUTSIDE
    if (dict % isPresent('transparent')) then
      call dict % get(matNames, 'transparent')
      allocate(mats(size(matNames) + 2))
      do i = 1, size(matNames)
        mats(i) = mm_nameMap % get(matNames(i))
      end do
      mats(i)   = OUTSIDE_MAT
      mats(i+1) = VOID_MAT
    else
      allocate(mats(2))
      mats(1) = OUTSIDE_MAT
      mats(2) = VOID_MAT
    end if
    
    ! Colourmap offset
    ! If not given select randomly
    if (dict % isPresent('offset')) then
      call dict % get(offset, 'offset')

    else
      call date_and_time(date, time)
      call FNV_1(date // time, offset)

    end if

    ! Img contains luminosity values, matIDs identifies which materials were hit
    call self % geom % rayPlot(lum, matIDs, camera, light, M, mats, fov, diffuse)

    ! Translate to an image.
    ! Obtain material colours and scale by luminosity
    matIDs = materialColour(matIDs, offset)
    img = matIDs

    ! Scale image by brightness
    img = brightnessScale(img, lum)

    ! Print image
    call imgBmp_toFile(img, outputFile)

  end subroutine makeRayPlot

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
  !! Scales an integer image by the provided value of brightness
  !!
  elemental function brightnessScale(img0, brightness) result(img)
    integer(shortInt), intent(in) :: img0
    real(defReal), intent(in)     :: brightness
    integer(shortInt)             :: img
    integer(shortInt)             :: r, g, b

    ! Extract RGB components (assuming 0x00RRGGBB)
    r = ishft(iand(img0, Z'00FF0000'), -16)
    g = ishft(iand(img0, Z'0000FF00'), -8)
    b = iand(img0, Z'000000FF')

    ! Scale and clamp to 0–255
    r = max(0, min(255, int(r * brightness)))
    g = max(0, min(255, int(g * brightness)))
    b = max(0, min(255, int(b * brightness)))

    ! Reassemble RGB integer
    img = ior(ior(ishft(r, 16), ishft(g, 8)), b)

  end function brightnessScale

end module visualiser_class
