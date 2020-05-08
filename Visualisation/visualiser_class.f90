!!
!! The visualiser object
!!
module visualiser_class

  use numPrecision
  use universalVariables
  use genericProcedures
  use commandLineUI,     only: getInputFile
  use dictionary_class
  use geometry_inter,    only: geometry
  use outputVTK_class

  implicit none
  private

  !!
  !! Object responsible for controlling visualisation
  !! Should be extensible for adding different visualisation methods
  !! Recieves and generates data for visualisation
  !!
  !! Requires a dictionary input which specifies the procedures to call
  !! Presently supports: VTK voxel mesh creation
  !!
  type, public :: visualiser
    character(nameLen)       :: name
    class(geometry), pointer :: geom => null()
  contains
    procedure :: init
    procedure :: makeVTK
  end type

contains

  !!
  !! Initialises visualiser and creates start-time
  !! visualisations
  !!
  subroutine init(self, geom, vizDict)
    class(visualiser), intent(inout)             :: self
    class(geometry), pointer, intent(inout)      :: geom
    class(dictionary), intent(in)                :: vizDict
    class(dictionary), pointer                   :: tempDict
    character(nameLen),dimension(:), allocatable :: keysArr
    integer(shortInt)                            :: i
    character(nameLen)                           :: type
    character(:), allocatable                    :: string
    character(nameLen) :: here ='init, visualiser_class.f90'

    ! Obtain file name
    call getInputFile(string)
    self % name = string

    ! Point to geometry
    self % geom => geom

    ! Loop through each dictionary and generate visualisation
    ! (if the visualisation method is available)
    call vizDict % keys(keysArr,'dict')

    do i=1,size(keysArr)
      tempDict => vizDict % getDictPtr(keysArr(i))
      call tempDict % get(type,'type')
      select case(type)
        case('vtk')
          call self % makeVTK(tempDict)
        case default
          call fatalError(here, 'Unrecognised visualisation - presently only accept vtk')
      end select

    end do

  end subroutine init

  !!
  !! Generate a VTK output
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
    call self % geom % voxelPlot(voxelMat, what, center, width)

    ! In principle, can add multiple data sets to VTK - not done here yet
    ! VTK data set will use 'what' variable as a name
    call vtk % addData(voxelMat, what)
    call vtk % output(self % name)
    call vtk % kill()

  end subroutine makeVTK

end module visualiser_class
