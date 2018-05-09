! Class which will contain all information regarding the geometry
! Contains:
!           array of surfaces
!           array of universes
!           array of cells
!           array of lattices
! Subroutines:
!           find whether neutron is in geometry & which cell it occupies
!
module geometry_class
  use numPrecision
  use genericProcedures
  use universalVariables

  use rng_class

  use dictionary_class
  use IOdictionary_class

  use coord_class

  use surface_class
  use sphere_class
  use plane_class
  use box_class
  use cylinder_class
  use squareCylinder_class
  use truncatedCylinder_class
  use cell_class
  use universe_class
  use lattice_class

  implicit none
  private

  type, public :: geometry
    type(surface_ptr), dimension(:), allocatable :: surfaces        ! pointers to all surfaces
    type(universe), dimension(:), allocatable    :: universes       ! array of all universes
    type(cell), dimension(:), allocatable        :: cells           ! array of all cells
    type(lattice), dimension(:), allocatable     :: lattices        ! array of all lattices
    type(surface_ptr)                            :: boundarySurface ! pointer to the boundary surface
    type(box)                                    :: boundingBox     ! bounding box for volume calculations
    type(universe_ptr)                           :: rootUniverse
    integer(shortInt)                            :: numSurfaces = 0
    integer(shortInt)                            :: numCells = 0
    integer(shortInt)                            :: numUniverses = 0
    integer(shortInt)                            :: numLattices = 0
  contains
    procedure :: init                    ! initialise geometry
    procedure :: whichCell               ! find which cell a neutron occupies
    procedure :: constructBoundingBox    ! Construct the box bounding the geometry
    procedure :: calculateVolumes        ! Calculates the volumes of all cells in the geometry
    procedure :: slicePlot               ! Produces a geometry plot of a slice perpendicular to a given axis
  end type geometry

contains

  subroutine init(self, surfPointers, cellArray, universeArray, rootIdx, boundarySurface, latticeArray)
    class(geometry), intent(inout)                     :: self
    type(surface_ptr), dimension(:)                    :: surfPointers
    class(cell),intent(in), dimension(:)               :: cellArray
    class(universe), intent(in), dimension(:)          :: universeArray
    class(lattice), intent(in), dimension(:), optional :: latticeArray
    integer(shortInt), intent(in)                      :: rootIdx
    class(surface_ptr)                                 :: boundarySurface

    allocate(self % surfaces (size(surfPointers)))
    allocate(self % cells (size(cellArray)))
    allocate(self % universes (size(universeArray)))

    self % surfaces = surfPointers
    self % cells = cellArray
    self % universes = universeArray
    self % numSurfaces = size(surfPointers)
    self % numCells = size(cellArray)
    self % numUniverses = size(universeArray)
    self % rootUniverse = self % universes(rootIdx)
    self % boundarySurface = boundarySurface

    if(present(latticeArray))then
      allocate(self % lattices (size(latticeArray)))
      self % lattices = latticeArray
      self % numLattices = size(latticeArray)
    end if

  end subroutine init

  !!
  !! Search universes until a base cell filled with material is located
  !!
  function whichCell(self, coords, n0) result(c)
    class(geometry), intent(in)             :: self
    class(coordList), intent(inout)         :: coords
    integer(shortInt), intent(in), optional :: n0
    real(defReal), dimension(3)             :: r, u, offset
    integer(shortInt), dimension(3)         :: ijkLat
    integer(shortInt)                       :: uniIdx, latIdx, n, instance
    integer(shortInt), dimension(3)         :: ijkUni
    type(cell_ptr)                          :: c
    type(universe_ptr)                      :: uni
    type(lattice_ptr)                       :: lat

    ! If a nesting level is not specified, begin from the root universe
    if(present(n0)) then
      n = n0
      uni = self % universes(coords % lvl(n) % uniIdx)
    else
      n = 1
      uni = self % rootUniverse
      call coords % resetNesting()
      coords % lvl(1) % uniIdx = uni % geometryIdx()
    end if

    ! Search through nested universe structures until a base cell is found
    do
      ! Set local co-ordinates
      r = coords % lvl(n) % r
      u = coords % lvl(n) % dir
      c = uni % whichCell(r, u)
      coords % lvl(n) % cellIdx = c % geometryIdx()

      ! If c is a null pointer then the point is not in the universe in which it was
      ! presumed to exist - this probably shouldn't happen!
      if (.NOT.c % associated()) then
        call fatalError('whichCell, geometry','The point could not be found')

      ! The cell is outside the geometry
      else if (.NOT. c % insideGeom()) then
        call uni % kill()
        call lat % kill()
        return

      ! The cell is filled with a universe - update the universe and continue search
      else if (c % fillType() == universeFill) then
        uniIdx = c % uniIdx()
        uni = self % universes(uniIdx)
        offset = uni % offset()
        call coords % addLevel(offset, uniIdx)
        n = n + 1
        cycle

      ! The cell is filled with a lattice - identify the lattice, move to the lattice co-ordinate
      ! system, then move to the lattice cell co-ordinate system
      else if (c % fillType() == latticeFill) then
        latIdx = c % latIdx()
        lat = self % lattices(latIdx)
        ijkLat = lat % findUniverse(r,u)
        uni = lat % universes(ijkLat)
        uniIdx = uni % geometryIdx()
        offset = uni % offset() + lat % localCoords(ijkLat)
        ijkIdx = lat % getijkIdx(ijkLat)
        call coords % addLevel(offset, uniIdx, latIdx, ijkIdx)
        n = n + 1
        cycle

      ! The cell is filled with a material - terminate the search successfully
      else if (c % fillType() == materialFill) then
        instance = c % coordCompare(coords)
        coords % matIdx = c % matIdx(instance)
        coords % uniqueCellID = c % uniqueID(instance)
        call uni % kill()
        call lat % kill()
        return

      ! Something has gone wrong...
      else
        print *, c % name()
        print *, c % materialIdx()
        print *, c % latIdx()
        print *, c % uniIdx()
        print *, c % geometryIdx()
        print *, c % fillType()
        call fatalError('whichCell, geometry','Could not find the cell')
      end if
    end do

  end function whichCell

  !!
  !! Assumes that the geometry is bounded by a cuboid
  !! Details of boundary box provided by user
  !!
  subroutine constructBoundingBox(self, orig, a)
    class(geometry), intent(inout)          :: self
    real(defReal), intent(in), dimension(3) :: orig
    real(defReal), intent(in), dimension(3) :: a
    call self % boundingBox % init(orig, a)
  end subroutine constructBoundingBox

  !!
  !! Calculate cell volumes by sampling numPoints points in the geometry
  !! subject to the geometry having a bounding box
  !!
  !! Will eventually require parallelisation
  !!
  subroutine calculateVolumes(self, numPoints, seed)
    class(geometry), intent(inout)           :: self
    integer(shortInt), intent(in)            :: numPoints
    integer(longInt), intent(in), optional   :: seed
    integer(longInt)                         :: s
    real(defReal), dimension(3)              :: minPoint, width, testPoint, dummyU
    real(defReal), dimension(:), allocatable :: cellVols
    integer(shortInt)                        :: i, cellIdx
    type(rng)                                :: rand
    type(box)                                :: bb
    type(cell_ptr)                           :: c
    type(coordList)                          :: coords

    allocate(cellVols(self % numCells))
    cellVols(:) = 0

    ! Cell searching requires a direction
    dummyU = [ONE, ZERO, ZERO]

    ! Replace by taking from clock?
    if(.not.present(seed)) then
      s = 1_8
    else
      s = seed
    end if
    call rand % init(s)

    bb = self % boundingBox
    minPoint = bb % origin - bb % a
    width = TWO * bb % a

    ! Volume calculation by point sampling
    do i = 1,numPoints
      testPoint(1) = minPoint(1) + width(1) * rand%get()
      testPoint(2) = minPoint(2) + width(2) * rand%get()
      testPoint(3) = minPoint(3) + width(3) * rand%get()

      call coords % init(testPoint, dummyU)

      c = self % whichCell(coords)
      cellIdx = c % geometryIdx()

      if (c % insideGeom()) then
        cellVols(cellIdx) = cellVols(cellIdx) + 1
      end if
    end do
    call c % kill()

    cellVols(:) = ONE*cellVols(:) * width(1) * width(2) * width(3) / numPoints
    print *,cellVols(:)
    ! Assign volumes to cells, dividing by the number of instances
    ! Only assign to base cells
    do i = 1,self%numCells
      if (self % cells(i) % fillType == materialFill) then
        self % cells(i) % volume = cellVols(i) / self % cells(i) % instances
      end if
    end do

  end subroutine calculateVolumes

  !!
  !! Routine to produce a geometry plot
  !! Given pixel specification, samples points uniformly in each direction
  !! Unique materals are assigned different colours
  !! Outputs a matrix of colours
  !! Specify:
  !! the number of pixels in each direction of the plane
  !! the point on which the plane is centred
  !! the (full) width in each direction of the plane
  !! the index of the perpendicular direction (optional - otherwise slice is XY)
  !!
  function slicePlot(self,pixels,centre,width,perpDirection)result(colourMatrix)
    class(geometry), intent(in)                    :: self
    integer(shortInt), dimension(2), intent(in)    :: pixels
    real(defReal), dimension(3), intent(in)        :: centre
    real(defReal), dimension(2), intent(in)        :: width
    integer(shortInt), dimension(:,:), allocatable :: colourMatrix
    integer(shortInt), optional, intent(in)        :: perpDirection
    integer(shortInt)                              :: perpD
    integer(shortInt)                              :: i,j
    real(defReal), dimension(3)                    :: x0, point, u
    type(cell_ptr)                                 :: c
    real(defReal), dimension(2)                    :: step
    type(coordList)                                :: coords

    allocate(colourMatrix(pixels(1),pixels(2)))

    ! Plot an x-y slice by default
    if(present(perpDirection)) then
      perpD = perpDirection
    else
      perpD = 3
    end if

    ! Create the co-ordinate points and widths
    step(:) = width(:) / pixels(:)

    ! Set a travel direction in case of surface intersections
    ! Also set the corner point x0
    if(perpD==1)then
      u = [ZERO, SQRT2_2, SQRT2_2]
      x0 = [centre(1), centre(2) - width(1)/TWO + step(1)/TWO, centre(3) - width(2)/TWO + step(2)/TWO]
    else if(perpD==2)then
      u = [SQRT2_2, ZERO, SQRT2_2]
      x0 = [centre(1) - width(1)/TWO + step(1)/TWO, centre(2), centre(3) - width(2)/TWO + step(2)/TWO]
    else if(perpD==3)then
      u = [SQRT2_2, SQRT2_2, ZERO]
      x0 = [centre(1) - width(1)/TWO + step(1)/TWO, centre(2) - width(2)/TWO + step(2)/TWO, centre(3)]
    else
      call fatalError('slicePlot, geometry','Incorrect plane direction supplied: must be between 1 and 3')
    end if

    ! Loop through pixels, finding base cells and identifying unique materials
    do i = 1,pixels(1)
      do j = 1,pixels(2)
        if (perpD==1) then
          point = [x0(1), x0(2) + i*step(1), x0(3) + j*step(2)]
        else if (perpD==2) then
          point = [x0(1) + i*step(1), x0(2), x0(3) + j*step(2)]
        else if (perpD==3) then
          point = [x0(1) + i*step(1), x0(2) + j*step(2), x0(3)]
        end if

        call coords % init(point, u)
        c = self % whichCell(coords)

        if (.not. c % insideGeom()) then
          colourMatrix(i,j) = 0
        else
          colourMatrix(i,j) = c % materialIdx()
        end if
      end do
    end do

    call c % kill()

  end function slicePlot

end module geometry_class
