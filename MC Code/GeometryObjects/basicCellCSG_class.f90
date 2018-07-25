module basicCellCSG_class

  use numPrecision
  use universalVariables
  use genericProcedures, only : fatalError
  use dictionary_class,  only : dictionary
  use nuclearData_inter, only : nuclearData

  ! Geometry modules
  !* Provisional
  use geometry_inter,    only : geometry
  use coord_class,       only : coord, coordList
  use csg_class,         only : csg, fillArray

  implicit none
  private

  type, public :: basicCellCSG
    private
    type(csg)  :: geom
  contains
    procedure :: init
    procedure :: placeCoord
    procedure :: whatIsAt
    procedure :: slicePlot

  end type basicCellCSG

contains

  !!
  !! Build basicCellCSG from dictionary and map of material names to matIdx
  !!
  subroutine init(self,dict, materials)
    class(basicCellCSG), intent(inout) :: self
    class(dictionary), intent(inout)   :: dict
    class(nuclearData), intent(in)     :: materials

    call self % geom % init(dict, materials)

  end subroutine init

  !!
  !! Places coordinate list into geometry
  !! Finds unique cell and material as well as coordinates at all intermediate levels
  !! Returns an error if coordList is uninitialised
  !!
  subroutine placeCoord(self, coords)
    class(basicCellCSG), intent(in) :: self
    type(coordList), intent(inout)  :: coords
    integer(shortInt)               :: nextUniRoot, fill
    integer(shortInt)               :: i, uniIdx
    real(defReal), dimension(3)     :: offset
    character(100), parameter  :: Here = 'placeCoord (basicCellCSG_class.f90)'


    ! Check that coordList is initialised
    if( coords % nesting < 1) call fatalError(Here,'Uninitialised coord List')

    associate (geom => self % geom )
      ! Place coordinates above geometry
      call coords % takeAboveGeom()

      ! Enter root universe
      coords % lvl(1) % uniIdx    = 1
      coords % lvl(1) % uniRootID = 1
      call geom % uShelf % shelf(1) % enter( coords % lvl(1), geom % cShelf, geom % sShelf)

      do i=1,hardcoded_max_nest
        ! Find cell fill
        call geom % fills % getFill(coords % lvl(i),fill, nextUniRoot)

        if (fill >= 0) then ! fill is a material or outside
          coords % matIdx = fill
          return

        end if

        if (fill < 0) then ! fill is a nested universe
          ! Change fill to uniIdx
          fill = -fill

          ! Current universe
          uniIdx = coords % lvl(i) % uniIdx

          ! Get cell offset
          offset = geom % uShelf % shelf(uniIdx) % cellOffset( coords % lvl(i) )

          ! Enter lower level
          call coords % addLevel(offset, fill, nextUniRoot)
          call geom % uShelf % shelf(fill)% enter( coords % lvl(i+1), geom % cShelf, geom % sShelf)

        end if
      end do

      call fatalError(Here,'Material cell was not found')

    end associate
  end subroutine placeCoord

  !!
  !! Given position in a geometry return material index and unique cell ID under r
  !!
  subroutine whatIsAt(self, r, matIdx, uniqueID)
    class(basicCellCSG), intent(in)        :: self
    real(defReal),dimension(3), intent(in) :: r
    integer(shortInt), intent(out)         :: matIdx
    integer(shortInt), intent(out)         :: uniqueId
    type(coordList)                        :: coords

    ! Initialise coordinates with +ve x direction
    call coords % init(r, [ONE, ZERO, ZERO])

    ! Place coordinates in geometry
    call self % placeCoord(coords)

    ! Read matIdx and uniqueId
    matIdx   = coords % matIdx
    uniqueId = coords % uniqueID()
    
  end subroutine whatIsAT

  !!
  !! Produce a 2D plot of a geometry
  !! Resolution is determined by a size of input matrix colorMat
  !! By default extent of a plot is determined by bounds of the domain and offset is [0,0,0]
  !!
  !! NOTES:
  !! -> what = {"material","cellID"} determines if matIdx is put in colorMat or unique cellID
  !! -> dir  = {"x","y","z"} specifies direction of the plot
  !! -> centre allows to set offset of a plane
  !! -> width sets well... width in each direction of the plane width(1) is either x or y
  !!
  subroutine slicePlot(self, colorMat, centre, dir, what, width)
    class(basicCellCSG), intent(in)                 :: self
    integer(shortInt),dimension(:,:), intent(inout) :: colorMat
    real(defReal), dimension(3), intent(in)         :: centre
    character(1), intent(in)                        :: dir
    character(*), intent(in)                        :: what
    real(defReal), dimension(2), optional           :: width
    real(defReal)                                   :: step1, step2
    real(defReal),dimension(3)                      :: x0, inc1, inc2, point
    integer(shortInt)                               :: i,j, flag, matIdx, uniqueId
    character(100),parameter :: Here = 'slicePlot (basicCellCSG_Class.f90)'

    ! Check that width was provided
    if (.not.present(width)) then
      call fatalError(Here,'Sorry but width must be provided before surfaces return their bounds')
    end if

    ! Calculate step in direction 1 & 2
    step1 = width(1) / size(colorMat,1)
    step2 = width(2) / size(colorMat,2)

    ! Depending on perpendicular direction set increment direction 1 & 2
    inc1 = ZERO
    inc2 = ZERO
    select case(dir)
      case('x')
        inc1(2) = ONE  ! Set inc1 to y-axis
        inc2(3) = ONE  ! Set inc2 to z-axis

      case('y')
        inc1(1) = ONE  ! Set inc1 to x-axis
        inc2(3) = ONE  ! Set inc2 to z-axis

      case('z')
        inc1(1) = ONE  ! Set inc1 to x-axis
        inc2(2) = ONE  ! Set inc2 to y-axis

      case default
        call fatalError(Here,'Unrecognised perpendicular ridection code: '//dir)
    end select

    ! Calculate corner
    x0 = centre - inc1 * (width(1)-step1)/2 - inc2 * (width(2)-step2)/2

    ! Choose to print uniqueID or matIdx
    select case(what)
      case('material')
       flag = 1

      case('uniqueID')
       flag = 2

      case default
        call fatalError(Here,'Target of plot must be material or uniqueID')
    end select

    ! Paint the color matrix (loop over leftmost index first for better memory efficiency
    do j = 1, size(colorMat,2)
      point = x0 + j * step2 * inc2

      do i = 1, size(colorMat,1)
        point = point + step1 * inc1

        ! Find material and unique id at the point
        call self % whatIsAt(point, matIdx, uniqueID)

        ! Paint the pixel of colorMat
        colorMat(i,j) = matIdx
        if (flag == 2) colorMat(i,j) = uniqueID

      end do
    end do
  end subroutine slicePlot


end module basicCellCSG_class
