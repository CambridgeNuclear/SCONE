module basicCellCSG_class

  use numPrecision
  use universalVariables
  use genericProcedures,  only : fatalError
  use vector_class,       only : vector
  use dictionary_class,   only : dictionary

  ! Map of material names to material indexes
  ! *** MAY be replaced with a char-int map later
  use nuclearData_inter,  only : nuclearData

  ! Geometry modules
  use cellGeometry_inter, only : cellGeometry
  use coord_class,        only : coord, coordList
  use csg_class,          only : csg, fillArray

  implicit none
  private

  !! Parameter for tracking recovery options
  integer(shortInt), parameter :: RECOVER = 0


  !!
  !! Implementation of the cell geometry model using CSG representation of geometry
  !! Does not support trip surface at the moment
  !!
  type, public,extends(cellGeometry) :: basicCellCSG
    private
    type(csg)  :: geom
  contains
    ! Build procedures
    procedure :: init

    ! Interface of geometry_inter
    procedure :: placeCoord
    procedure :: whatIsAt
    procedure :: bounds
    procedure :: slicePlot
    procedure :: voxelPlot

    ! Interface of cellGeometry_inter
    procedure :: move
    procedure :: teleport
    procedure :: moveGlobal

    ! Private procedures
    procedure, private :: closestDist
    procedure, private :: diveToMat

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

      ! Dive into geometry until material is found from 1st level
      call self % diveToMat( coords, 1)

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
  !! Return bounds of the geometry domain in cartesian co-ordinates
  !! bounds = [x_min, x_max, y_min, y_max, z_min, z_max]
  !! If geometry is infinate in a given axis direction * then:
  !! *_min = *_max = ZERO
  !!
  function bounds(self)
    class(basicCellCSG), intent(in) :: self
    real(defReal),dimension(6)      :: bounds
    character(100), parameter :: Here ='bounds (basicCellCSG_class.f90)'

    call fatalError(Here,'Bounds are not yet implemented. Need to modify surfaces to provide them')

  end function bounds

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
    class(basicCellCSG), intent(in)                  :: self
    integer(shortInt),dimension(:,:), intent(inout)  :: colorMat
    real(defReal), dimension(3), intent(in)          :: centre
    character(1), intent(in)                         :: dir
    character(*), intent(in)                         :: what
    real(defReal), dimension(2), optional,intent(in) :: width
    real(defReal)                                    :: step1, step2
    real(defReal),dimension(3)                       :: x0, inc1, inc2, point
    integer(shortInt)                                :: i,j, flag, matIdx, uniqueId
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
    x0 = centre - inc1 * (width(1) + step1)/2 - inc2 * (width(2) + step2)/2

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

  !!
  !! Produce a voxel 3D plot of geometry
  !! Resolution is determined by a size of input voxelMat
  !! By default bounds of the voxel plot correspond to bounds of geometry
  !!
  !! NOTES:
  !! -> what = {"material","cellID"} determines if matIdx is put in voxel Mat or unique cellID
  !! -> centere and optional width specify extent of plot
  !! -> Voxel plot is always a box and it is axis aligned
  !!
  subroutine voxelPlot(self,voxelMat,what,center,width)
    class(basicCellCSG),intent(in)                   :: self
    integer(shortInt),dimension(:,:,:),intent(inout) :: voxelMat
    character(*),intent(in)                          :: what
    real(defReal),dimension(3),intent(in)            :: center
    real(defReal),dimension(3),optional,intent(in)   :: width
    real(defReal)                                    :: stepZ, z0
    real(defReal),dimension(3)                       :: point
    real(defReal),dimension(2)                       :: width_bar
    integer(shortInt)                                :: i
    character(100),parameter :: Here = 'voxelPlot ( basicCellCSG_class.f90)'

    ! Check that width was provided
    if (.not.present(width)) then
      call fatalError(Here,'Sorry but width must be provided before surfaces return their bounds')
    end if

    ! Calculate step in z direction
    stepZ = width(3) / size(voxelMat,3)

    ! Calculate starting z
    z0 = center(3) - (width(3) + stepZ)/TWO

    ! Construct voxel plot from multiple slice plots
    point = center
    do i=1,size(voxelMat,3)
      point(3) = z0 + stepZ
      call self % slicePlot(voxelMat(:,:,i), point, 'z', what, width(1:2) )

    end do

 end subroutine voxelPlot

  !!
  !! Given coordinates placed in a geometry move point through the geometry
  !! Move by maxDist. Stop at the boundaries between unique cell IDs or at the trip surface
  !!
  !! If during motion boundary is crossed transformations are applied and particle is stopped
  !! If particle is stopped at the trip surface tripFlag is set to .true.
  !! If particle stoped at the collision isColl is set to .true.

  !! NOTE:
  !!  -> if coords is not placed in the geometry behaviour is unspecified
  !!  -> if maxDist < 0.0 behaviour is unspecified
  !!
  subroutine move(self, coords, maxDist, isColl, tripFlag)
    class(basicCellCSG), intent(inout)    :: self
    type(coordList), intent(inout)        :: coords
    real(defReal),intent(in)              :: maxDist
    logical(defBool),intent(out)          :: isColl
    logical(defBool),optional,intent(out) :: tripFlag
    integer(shortInt)                     :: surfIdx, level, uniIdx
    real(defReal)                         :: dist
    type(vector)                          :: r, u
    character(100),parameter  :: Here ='move (basicCellCSG_Class.f90)'

    ! *** TRIP SURFACES NOT YET IMPLEMENTED
    if(present(tripFlag)) tripFlag = .false.

    ! Perform Check whether coordinates are placed
    if(.not. coords % isPlaced()) then
      call fatalError(Here, 'Coordinate List is not placed in geometry')
    end if

    ! Find distance to the next surface
    call self % closestDist(dist, surfIdx, level, coords)

    !print *, surfIdx, dist

    if (maxDist < dist) then ! Moves within cell
      call coords % moveLocal(maxDist, coords % nesting)
      isColl = .true.

    else if (surfIdx == self % geom % boundaryIdx) then ! Crosses domain boundary
      ! Move global to the boundary
      call coords % moveGlobal(dist)

      ! Apply boundary condition *** This interface is to horrible. Must be repaired
      r = coords % lvl(1) % r
      u = coords % lvl(1) % dir

      call self % geom % sShelf % shelf(surfIdx) % boundaryTransform(r, u)

      coords % lvl(1) % r   = r % v
      coords % lvl(1) % dir = u % v

      ! Return particle to geometry
      call self % placeCoord(coords)
      isColl = .false.

    else ! Crosses cell to cell boundary
      ! Move to the boundary at "level"
      call coords % moveLocal(dist, level)
      uniIdx = coords % lvl(level) % uniIdx
      isColl = .false.

      associate (g => self % geom)

        ! Cross boundary (sets cellIdx and localCell after boundary)
        call g % uShelf % shelf(uniIdx) % cross(coords % lvl(level), surfIdx, g % cShelf, g % sShelf)

        ! Dive further into geometry until material is found
        call self % diveToMat(coords,level)

      end associate

    end if
  end subroutine move

  !!
  !! Given coordinates transport point outside the geometry by distance dist
  !!
  !! If during motion boundary is crossed transformation are applied
  !! Particle is NOT stoped but transport continues until total
  !! transport length is equal to dist.
  !!
  !! NOTE:
  !!  -> if coords is not placed in the geometry behaviour is unspecified
  !!  -> if maxDist < 0.0 behaviour is unspecified
  !!
  subroutine teleport(self,coords,dist)
    class(basicCellCSG), intent(inout) :: self
    type(coordList), intent(inout)     :: coords
    real(defReal), intent(in)          :: dist
    type(vector)                       :: r, u

    ! Move the coords above the geometry
    call coords % moveGlobal(dist)

    ! Place coordinates in back in the goemetry
    call self % placeCoord(coords)

    ! If point is outside apply boundary transformations
    if (coords % matIdx == OUTSIDE_FILL) then

      ! Apply boundary condition *** This interface is to horrible. Must be repaired
      r = coords % lvl(1) % r
      u = coords % lvl(1) % dir

      call self % geom % sShelf % shelf(self % geom % boundaryIdx) % boundaryTransform(r, u)

      coords % lvl(1) % r   = r % v
      coords % lvl(1) % dir = u % v

      ! Return particle to geometry
      call self % placeCoord(coords)

    end if

  end subroutine teleport

  !!
  !! Given coordinates move point in global geometry level
  !! Move by maxDist. Stop at the boundary or at the trip surface
  !!
  !! If during motion boundary is crossed transformations are applied and particle is stopped
  !! If particle is stopped at the trip surface tripFlag is set to .true.
  !!
  !! NOTE:
  !!  -> if maxDist < 0.0 behaviour is unspecified
  !!
  subroutine moveGlobal(self,coords,maxDist,tripFlag)
    class(basicCellCSG), intent(inout)    :: self
    type(coordList), intent(inout)        :: coords
    real(defReal), intent(in)             :: maxDist
    logical(defBool),optional,intent(out) :: tripFlag
    integer(shortInt)                     :: surfIdx
    real(defReal)                         :: dist
    type(vector)                          :: r, u
    character(100), parameter :: Here ='moveGlobal (basciCSG_class.f90)'

    ! *** TRIP SURFACES NOT YET IMPLEMENTED
    if(present(tripFlag)) tripFlag = .false.

    ! Perform Check whether coordinates are placed
    if(.not. coords % isPlaced()) then
      call fatalError(Here, 'Coordinate List is not placed in geometry')
    end if

    ! Find distance to the boundary surface
    associate ( g => self % geom)
      surfIdx = g % boundaryIdx

      ! Get position and direction into vectors
      r = coords % lvl(1) % r
      u = coords % lvl(1) % dir

      ! Calculate distance to the boundary
      call g % sShelf % shelf(surfIdx) % distance(dist, surfIdx, r, u)

      if (maxDist < dist) then
        ! Move above the geometry
        call coords % moveGlobal(maxDist)

      else
        ! Move above the geometry
        call coords % moveGlobal(dist)

        ! Apply boundary condition *** This interface is to horrible. Must be repaired
        r = coords % lvl(1) % r
        u = coords % lvl(1) % dir

        call g % sShelf % shelf(g % boundaryIdx) % boundaryTransform(r, u)

        coords % lvl(1) % r   = r % v
        coords % lvl(1) % dir = u % v
      end if

      ! Return particle to geometry
      call self % placeCoord(coords)

    end associate
  end subroutine moveGlobal


  !!
  !! Returns distance, surfIdx for the next surface
  !! Also returns level at which the surface is encountered
  !!
  subroutine closestDist(self, dist, surfIdx, level, coords)
    class(basicCellCSG), intent(in)       :: self
    real(defReal), intent(out)            :: dist
    integer(shortInt), intent(out)        :: surfIdx
    integer(shortInt), intent(out)        :: level
    type(coordList), intent(in)           :: coords
    integer(shortInt)                     :: l, sIdxTest, uniIdx
    real(defReal)                         :: testDist

    dist = INFINITY
    associate ( g => self % geom)
      ! Loop over levels occupied by the particle
      do l=1, coords % nesting
        ! Get universe index
        uniIdx = coords % lvl(l) % uniIdx

        ! Find distance to next surface at level l
        call g % uShelf % shelf(uniIdx) % distance(testDist       ,&
                                                   sIdxTest       ,&
                                                   coords % lvl(l),&
                                                   g % cShelf     ,&
                                                   g % sShelf     )

        ! Save distance, surfIdx & level coresponding to shortest distance
        if (testDist < dist) then
          dist    = testDist
          surfIdx = sIdxTest
          level   = l
        end if
      end do
    end associate
  end subroutine closestDist

  !!
  !! Dives into geometry until material is found starting from level start
  !! Assumes that coords are already placed at level start
  !!
  subroutine diveToMat(self, coords, start)
    class(basicCellCSG), intent(in) :: self
    type(coordList), intent(inout)  :: coords
    integer(shortInt), intent(in)   :: start
    integer(shortInt)               :: uniIdx, fill, nextUniRoot
    integer(shortInt)               :: i
    real(defReal),dimension(3)      :: offset
    character(100), parameter :: Here ='diveToMat (basicCellCSG_class.f90)'

    associate (g => self % geom)

      do i=start,hardcoded_max_nest
        ! Find cell fill
        call g % fills % getFill(coords % lvl(i),fill, nextUniRoot)

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
          offset = g % uShelf % shelf(uniIdx) % cellOffset( coords % lvl(i) )

          ! Enter lower level
          call coords % addLevel(offset, fill, nextUniRoot)
          call g % uShelf % shelf(fill)% enter( coords % lvl(i+1), g % cShelf, g % sShelf)

        end if
      end do

      call fatalError(Here,'Material cell was not found')

    end associate
  end subroutine diveToMat

end module basicCellCSG_class
