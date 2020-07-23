module geometry_inter

  use numPrecision
  use universalVariables, only : X_AXIS, Y_AXIS, Z_AXIS
  use genericProcedures,  only : fatalError
  use dictionary_class,   only : dictionary
  use charMap_class,      only : charMap
  use coord_class,        only : coordList

  implicit none
  private

  !!
  !! Abstract interface for all geometry implementation
  !!
  type, public, abstract :: geometry
  contains
    ! Deferred procedures
    procedure(init), deferred       :: init
    procedure(kill), deferred       :: kill
    procedure(placeCoord), deferred :: placeCoord
    procedure(whatIsAt), deferred   :: whatIsAt
    procedure(bounds), deferred     :: bounds
    procedure(move), deferred       :: move
    procedure(moveGlobal), deferred :: moveGlobal
    procedure(teleport), deferred   :: teleport

    ! Common procedures
    procedure :: slicePlot
    procedure :: voxelplot
  end type geometry

  abstract interface

    !!
    !! Initialise geometry
    !!
    !! Args:
    !!   dict [in] -> Dictionary with geometry definition
    !!   mats [in] -> Map of material names to matIdx
    !!   silent [in] -> Optional. Set to .true. to surpress console messeges. Default .false.
    !!
    subroutine init(self, dict, mats, silent)
      import :: geometry, dictionary, charMap, defBool
      class(geometry), intent(inout)         :: self
      class(dictionary), intent(in)          :: dict
      type(charMap), intent(in)              :: mats
      logical(defBool), optional, intent(in) :: silent
    end subroutine init

    !!
    !! Return to uninitialised state
    !!
    elemental subroutine kill(self)
      import :: geometry
      class(geometry), intent(inout) :: self
    end subroutine kill

    !!
    !! Place coordinate list into geometry
    !!
    !! Finds unique cell and material as well as location at all intermediate levels
    !!
    !! Args:
    !!   coords [inout] -> Initialised coordinate list. This means that location in tope level must
    !!     be valid and direction normalised to 1.0.
    !!
    !! Errors:
    !!   fatalError if coordList is not initialised
    !!
    subroutine placeCoord(self, coords)
      import :: geometry, coordList
      class(geometry), intent(in)    :: self
      type(coordList), intent(inout) :: coords
    end subroutine placeCoord

    !!
    !! Find material and unique cell at a given location
    !!
    !! Args:
    !!   matIdx [out] -> material index at the location
    !!   uniqueID [out] -> Unique Id at the location
    !!   r [in] -> Position in the geometry
    !!   u [in] -> Optional. Normalised direction (norm2(u) = 1.0) (default = [1, 0, 0])
    !!
    subroutine whatIsAt(self, matIdx, uniqueID, r, u)
      import :: geometry, shortInt, defReal
      class(geometry), intent(in)                       :: self
      integer(shortInt), intent(out)                    :: matIdx
      integer(shortInt), intent(out)                    :: uniqueID
      real(defReal), dimension(3), intent(in)           :: r
      real(defReal), dimension(3), optional, intent(in) :: u
    end subroutine whatIsAt

    !!
    !! Return Axis Aligned Bounding Box encompassing the geometry
    !!
    !! Provides with bounds of the geometry
    !!
    !! Args:
    !!   None
    !!
    !! Result:
    !!   Size 6 array [x_min, y_min, z_min, x_max, y_max, z_max] with locations of
    !!   the lower and the high corner of the axis aligned bounding box.
    !!   If geometry is infinate in a given axis direction * then *_min = *_max = ZERO
    !!
    function bounds(self)
      import :: geometry, defReal
      class(geometry), intent(in) :: self
      real(defReal), dimension(6) :: bounds
    end function bounds

    !!
    !! Given coordinates placed in the geometry move point through the geometry
    !!
    !! Move by up to maxDist stopping at domain boundary or untill matIdx or uniqueID changes
    !! When particle hits boundary, boundary conditions are applied before returning.
    !!
    !! Following events can be returned:
    !!   COLL_EV      -> Particle moved by entire maxDist. Collision happens
    !!   BOUNDARY_EV  -> Particle hit domain boundary
    !!   CROSS_EV     -> Partilce crossed to a region with different matIdx or uniqueID
    !!   LOST_EV      -> Something gone wrong in tracking and particle is lost
    !!
    !! Args:
    !!   coords [inout]  -> Coordinate list of the particle to be moved through the geometry
    !!   maxDict [inout] -> Maximum distance to move the position. If movment is stopped
    !!     prematurely (e.g. hitting boundary), maxDist is set to the distance the particle has
    !!     moved by.
    !!   event [out] -> Event flag that specifies what finished the movement.
    !!
    !! Errors:
    !!   If coords is not placed in the geometry behaviour is unspecified
    !!   If maxDist < 0.0 behaviour is unspecified
    !!
    subroutine move(self, coords, maxDist, event)
      import :: geometry, coordList, defReal, shortInt
      class(geometry), intent(in)    :: self
      type(coordList), intent(inout) :: coords
      real(defReal), intent(inout)   :: maxDist
      integer(shortInt), intent(out) :: event
    end subroutine move

    !!
    !! Move a particle in the top (global) level in the geometry
    !!
    !! Move up to maxDist or untill domain boundary is hit, in which case applies boundary
    !! conditions and exits.
    !!
    !! Following events can be returned:
    !!   COLL_EV      -> Particle moved by entire maxDist. Collision happens
    !!   BOUNDARY_EV  -> Particle hit domain boundary
    !!
    !! Args:
    !!   coords [inout] -> Initialised (but not necesserly placed) coordList for a particle to be
    !!     moved. Will become placed on exit.
    !!   maxDict [inout] -> Maximum distance to move the position. If movment is stopped
    !!     prematurely (e.g. hitting boundary), maxDist is set to the distance the particle has
    !!     moved by.
    !!   event [out] -> Event flag that specifies what finished the movement.
    !!
    !! Errors:
    !!   If maxDist < 0.0 behaviour is unspecified
    !!
    subroutine moveGlobal(self, coords, maxDist, event)
      import :: geometry, coordList, defReal, shortInt
      class(geometry), intent(in)    :: self
      type(coordList), intent(inout) :: coords
      real(defReal), intent(inout)   :: maxDist
      integer(shortInt), intent(out) :: event
    end subroutine moveGlobal

    !!
    !! Move a particle in the top level without stopping
    !!
    !! Moves exactly by a given distance. If domain boundary is hit, boundary conditions are
    !! applied and movement continious untill full distance is reached.
    !!
    !! Args:
    !!   coords [inout] -> Initialised (but not necesserly placed) coordList for a particle to be
    !!     moved. Will become placed on exit.
    !!   dist [in] -> Distance by which move the particle
    !!
    !! Errors:
    !!   If maxDist < 0.0 behaviour is unspecified
    !!
    subroutine teleport(self, coords, dist)
      import :: geometry, coordList, defReal
      class(geometry), intent(in)    :: self
      type(coordList), intent(inout) :: coords
      real(defReal), intent(in)      :: dist
    end subroutine teleport

  end interface

contains

  !!
  !! Produce a 2D plot of the geometry
  !!
  !! Resolution is determined by a size of provided output matrix
  !! By default plot plane is normal to z-axis, witch width determined by bounds of the
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
    class(geometry), intent(in)                       :: self
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
    character(100), parameter :: Here = 'slicePlot (geometry_inter.f90)'

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
      aabb = self % bounds()
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

    do j = 1, size(img, 2)
      point(plane(2)) = corner(plane(2)) + step(plane(2)) * j

      do i = 1, size(img, 1)
        point(plane(1)) = corner(plane(1)) + step(plane(1)) * i

        ! Find material and paint image
        call self % whatIsAt(matIdx, uniqueID, point)

        ! Paint the pixel
        if (printMat) then
          img(i, j) = matIdx
        else
          img(i, j) = uniqueID
        end if

      end do
    end do

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
    class(geometry), intent(in)                       :: self
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
      aabb = self % bounds()
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


end module geometry_inter
