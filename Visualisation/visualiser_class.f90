!!
!! The visualiser object  
!!
module visualiser_class

  use numPrecision
  use universalVariables
  use genericProcedures, only: fatalError
  use commandLineUI,     only: getInputFile
  use dictionary_class,  only: dictionary
  use geometry_inter,    only: geometry
  use outputVTK_class

  implicit none
  private

  !!
  !! Object responsible for controlling visualisation
  !!
  !! Object that creates images relating to SCONE geometries
  !! Should be extensible for adding different visualisation methods
  !! Recieves and generates data for visualisation
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
  !!   makeVTK -> constructs VTK files
  !!   kill    -> cleans up visualiser
  !!
  !! Sample dictionary input:
  !!   viz{
  !!     vizDict1{ <outputVTK> }
  !!     #vizDict2{ <outputVTK> }#
  !!   }
  !!
  type, public :: visualiser
    character(nameLen), private       :: name
    class(geometry), pointer, private :: geom => null()
    type(dictionary), private         :: vizDict
  contains
    procedure :: init
    procedure :: makeViz
    procedure :: makeVTK
    procedure :: kill
  end type

contains

  !!
  !! Initialises visualiser
  !!
  !! Provides visualiser with filename for output,
  !! geometry information, and the dictionary decribing
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
        case default
          call fatalError(here,&
           'Unrecognised visualisation - presently only accept vtk')
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
  !! Result:
  !!   A vtk visualisation
  !!
  !! Errors:
  !!   Returns an error if there is an incorrect size for any of the 
  !!   required vtk inputs
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

    if (size(corner) .NE. 3) then
      call fatalError(here,'Voxel plot requires corner to have 3 values')
    endif
    if (size(width) .NE. 3) then
      call fatalError(here,'Voxel plot requires width to have 3 values')
    endif
    if (size(nVox) .NE. 3) then
      call fatalError(here,'Voxel plot requires vox to have 3 values')
    endif
    allocate(voxelMat(nVox(1),nVox(2),nVox(3)))

    ! Have geometry obtain data
    call self % geom % voxelPlot(voxelMat,what,center,width)

    ! In principle, can add multiple data sets to VTK - not done here yet
    ! VTK data set will use 'what' variable as a name
    call vtk % addData(voxelMat,what)
    call vtk % output(self % name)
    call vtk % kill()

  end subroutine makeVTK

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

end module visualiser_class

