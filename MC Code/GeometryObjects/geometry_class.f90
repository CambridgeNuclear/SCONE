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
  use universalVariables
  use genericProcedures, only : fatalError


  use particle_class,       only : particle

  use RNG_class,          only : RNG

  use dictionary_class,   only : dictionary

  use nuclearData_inter,  only : nuclearData

  use coord_class,        only : coord, coordList

  use surface_inter,       only : surface_ptr
  use surfaceFactory_func, only : new_surface_ptr


  use box_class,           only : box

  use cell_class,          only : cell, cell_ptr
  use universe_class,      only : universe, universe_ptr
  use lattice_class,       only : lattice, lattice_ptr

  implicit none
  private

  type, public :: geometry
    type(surface_ptr), dimension(:), allocatable :: surfaces             ! pointers to all surfaces
    type(universe), dimension(:), allocatable    :: universes            ! array of all universes
    type(cell), dimension(:), allocatable        :: cells                ! array of all cells
    type(lattice), dimension(:), allocatable     :: lattices             ! array of all lattices
    type(surface_ptr)                            :: boundarySurface      ! pointer to the boundary surface
    type(box)                                    :: boundingBox          ! bounding box for volume calculations
    logical(defBool)                             :: boxDefined = .FALSE. ! logical to check if the bounding box has been defined
    type(universe_ptr)                           :: rootUniverse
    integer(shortInt)                            :: numSurfaces = 0
    integer(shortInt)                            :: numCells = 0
    integer(shortInt)                            :: numUniverses = 0
    integer(shortInt)                            :: numLattices = 0
    integer(shortInt)                            :: numRegions = 0
  contains
    ! Build Procedures
    procedure :: init                    ! initialise geometry
    procedure :: initSurfaces            ! initialise surfaces
    procedure :: initCells               ! initialise cells
    procedure :: initUniverses           ! initialise universes
    procedure :: initLattices            ! initialise lattices
    procedure :: fillCells               ! fill cells with non-material contents
    procedure :: cellFromDict            ! initialise a cell from a dictionary
    procedure :: universeFromDict        ! initialise a universe from a dictionary
    procedure :: latticeFromDict         ! initialise a lattice from a dictionary
    procedure :: fillCellUni             ! fill a cell with a universe
    procedure :: fillCellLat             ! fill a cell with a lattice
    procedure :: setBoundaryConditions   ! set the boundary conditions for the geometry
    procedure :: fillMaterialCells       ! fill cells with materials

    ! ? Procedures
    procedure :: ennumerateRegions       ! Find all unique regions and assign an ID
    procedure :: countCell               ! Recursively count all cell instances
    procedure :: locateCell              ! Recursively locate all cell instances
    procedure :: constructBoundingBox    ! Construct the box bounding the geometry

    ! Interface Procedures
    procedure :: whichCell               ! find which cell a neutron occupies
    procedure :: calculateVolumes        ! Calculates the volumes of all cells in the geometry
    procedure :: slicePlot               ! Produces a geometry plot of a slice perpendicular to a given axis
    procedure :: voxelPlot               ! Produces a voxel plot of a chosen section of the geometry
    procedure :: placeParticle           ! Places a fresh particle in geometry. Assigns region and mat ID

  end type geometry

contains


  !!
  !! Puts fresh particle in a geometry
  !! Assigns regionID and matIdx to a particle
  !!
  subroutine placeParticle(self,p)
    class(geometry), intent(in)   :: self
    type(particle), intent(inout) :: p
    type(cell_ptr)                :: dummy_c

    dummy_c = self % whichCell(p % coords)
    call p % updateLocation()

  end subroutine placeParticle

  !!
  !! Initialise geometry from a dictionary
  !!
  subroutine init(self, dict, materials)
    class(geometry), intent(inout)     :: self
    class(dictionary), intent(inout)   :: dict
    class(nuclearData), intent(inout)  :: materials
    type(dictionary)                   :: surfDict, cellDict, uniDict, latDict

    print *,'Reading geometry from dictionary'

    surfDict = dict % getDict('surfaces')
    call self % initSurfaces(surfDict)
    call surfDict % kill()

    cellDict = dict % getDict('cells')
    call self % initCells(cellDict)

    uniDict = dict % getDict('universes')
    call self % initUniverses(uniDict)
    call uniDict % kill()

    if(dict % isPresent('lattices')) then
      latDict = dict % getDict('lattices')
      call self % initLattices(latDict)
      call latDict % kill()
    else
      self % numLattices = 0
    end if

    print *,'Filling non-material cells'
    call self % fillCells(cellDict, dict)

    print *,'Identifying unique material cell instances'
    call self % ennumerateRegions()

    print *,'Filling material cells'
    call self % fillMaterialCells(cellDict, materials)

    call cellDict % kill()
    call dict % kill()

  end subroutine init

  !!
  !! Initialise the surfaces of which the geometry is composed
  !!
  subroutine initSurfaces(self, surfDict)
    class(geometry), intent(inout)                 :: self
    class(dictionary), intent(in)                  :: surfDict
    character(nameLen), dimension(:), allocatable  :: keys
    type(dictionary)                               :: tempDict
    integer(shortInt)                              :: i, j, id, testId

    call surfDict % keysDict(keys)
    self % numSurfaces = size(keys)
    allocate(self % surfaces(self % numSurfaces))

    ! Construct surfaces and store in the surface array
    print *,'Constructing ',self % numSurfaces,' surfaces'
    do i=1,self % numSurfaces
      call surfDict % get(tempDict,keys(i))
      self % surfaces(i) = new_surface_ptr(tempDict,keys(i))

    end do

    ! Ensure no surfaces have the same id's
    if (self % numSurfaces > 1) then
      do i=1,self % numSurfaces-1
        id = self % surfaces(i) % id()
        do j=i+1,self % numSurfaces
          testId = self % surfaces(j) % id()
          if(testId==id) call fatalError('initSurfaces, geometry','Surfaces defined with identical ids')
        end do
      end do
    end if

  end subroutine initSurfaces

  !!
  !! Initialise cells from dictionary
  !!
  subroutine initCells(self, cellDict)
    class(geometry), intent(inout)                :: self
    class(dictionary), intent(in)                 :: cellDict
    character(nameLen), dimension(:), allocatable :: keys
    integer(shortInt)                             :: outsideDefined, i, j, id, testId

    call cellDict % keysDict(keys)
    self % numCells = size(keys)
    allocate(self % cells(self % numCells))

    ! Construct cells and store in the cell array
    ! Also ensure that the outside cell is defined
    outsideDefined = 0
    print *,'Constructing ',self % numCells,' cells'
    do i=1,self % numCells
      call self % cellFromDict(cellDict % getDict(keys(i)), i, keys(i))
      if (self % cells(i) % fillType == outsideFill) then
        outsideDefined = outsideDefined + 1
      end if
    end do
    deallocate(keys)

    if (outsideDefined < 1) then
      call fatalError('initCells, geometry','Must define one outside cell')
    else if (outsideDefined > 1) then
      call fatalError('initCells, geometry','More than one outside cell defined')
    end if

    if (self % numCells < 2) then
      call fatalError('initCells, geometry','Must define at least an outside and a material cell')
    end if

    ! Ensure no cells have the same id's
    do i=1,self % numCells-1
      id = self % cells(i) % id
      do j=i+1,self % numCells
        testId = self % cells(j) % id
        if(testId==id) call fatalError('initCells, geometry','Cells defined with identical ids')
      end do
    end do

  end subroutine initCells

  !!
  !! Initialise universes from dictionary
  !!
  subroutine initUniverses(self, uniDict)
    class(geometry), intent(inout)                :: self
    class(dictionary), intent(in)                 :: uniDict
    character(nameLen), dimension(:), allocatable :: keys
    integer(shortInt)                             :: i, j, id, testId, rootIdx
    logical(defBool)                              :: foundRoot
    type(dictionary)                              :: uniDaughterDict

    ! Read how many universes exist
    call uniDict % keysDict(keys)
    self % numUniverses = size(keys)
    allocate(self % universes(self % numUniverses))

    ! Construct universes and store in the universe array
    print *,'Constructing ',self % numUniverses,' universes'
    foundRoot = .FALSE.
    do i=1,self % numUniverses
      uniDaughterDict = uniDict % getDict(keys(i))
      if(uniDaughterDict % isPresent('root')) then
        if (uniDaughterDict % getInt('root') == 1) then
          call self % universeFromDict(uniDaughterDict, i, keys(i), .TRUE.)
          rootIdx = i
          if (.NOT.foundRoot) then
            foundRoot = .TRUE.
          else
            call fatalError('initUniverses, geometry','More than one root universe defined')
          end if
        end if
      else
        call self % universeFromDict(uniDaughterDict, i, keys(i), .FALSE.)
      end if
    end do

    ! Ensure that a root universe has been defined. Store pointer to it.
    if(.not.foundRoot) call fatalError('initUniverses, geometry','The root universe has not been defined')
    self % rootUniverse = self % universes(rootIdx) !*** This assignment is not STANDARD FORTRAN. self % universes are not "target"

    ! Ensure no universes have the same id's
    if (self % numUniverses > 1) then
      do i=1,self % numUniverses-1
        id = self % universes(i) % id
        do j=i+1,self % numUniverses
          testId = self % universes(j) % id
          if(testId==id) call fatalError('initUniverses, geometry',&
          'Universes defined with identical ids')
        end do
      end do
    end if

  end subroutine initUniverses

  !!
  !! Initialise lattices from dictionary
  !!
  subroutine initLattices(self, latDict)
    class(geometry), intent(inout)                :: self
    class(dictionary), intent(in)                 :: latDict
    character(nameLen), dimension(:), allocatable :: keys
    integer(shortInt)                             :: i, j, id, testId

    call latDict % keysDict(keys)
    self % numLattices = size(keys)
    allocate(self % lattices(self % numLattices))

    print *,'Constructing ',self % numLattices,' lattices'

    ! Construct lattices and store in the lattice array
    do i=1,self % numLattices
      call self % latticeFromDict(latDict % getDict(keys(i)), keys(i), i)
    end do
    deallocate(keys)

    ! Ensure no lattices have the same id's
    if (self % numLattices > 1) then
      do i=1,self % numLattices-1
        id = self % lattices(i) % id
        do j=i+1,self % numLattices
          testId = self % lattices(j) % id
          if(testId==id) call fatalError('initLattices, geometry','Lattices defined with identical ids')
        end do
      end do
    end if

  end subroutine initLattices

  subroutine fillCells(self, cellDict, geomDict)
    class(geometry), intent(inout)                :: self
    class(dictionary), intent(in)                 :: cellDict
    class(dictionary), intent(in)                 :: geomDict
    character(nameLen), dimension(:), allocatable :: keys
    integer(shortInt)                             :: i, j, fillType

    call cellDict % keysDict(keys)
    print *,'Filling cells with universes and lattices'
    do i=1,self % numCells
      fillType = self % cells(i) % fillType
      if (fillType == universeFill) then
        call self % fillCellUni(i, cellDict % getDict(keys(i)))
      else if (fillType == latticeFill) then
        call self % fillCellLat(i, cellDict % getDict(keys(i)))
      else if (fillType == materialFill) then
        cycle ! Fill with materials after number of cell instances have been identified
      else if (fillType == outsideFill) then
        ! Read and apply boundary conditions
        print *,'Applying boundary conditions'
        call self % setBoundaryConditions(i, geomDict)
        cycle
      else
        call fatalError('fillCell, geometry','Cell has an incorrect fill type')
      end if
    end do

  end subroutine fillCells

  !!
  !! Search universes until a base cell filled with material is located
  !!
  function whichCell(self, coords, n0) result(c)
    class(geometry), intent(in)             :: self
    class(coordList), intent(inout)         :: coords
    integer(shortInt), intent(in), optional :: n0
    real(defReal), dimension(3)             :: r, u, offset
    integer(shortInt), dimension(3)         :: ijkLat
    integer(shortInt)                       :: uniIdx, latIdx, ijkIdx, n, instance
    integer(shortInt), dimension(3)         :: ijkUni
    type(cell_ptr)                          :: c
    type(universe_ptr)                      :: uni
    type(lattice_ptr)                       :: lat

    ! If a nesting level is not specified, begin from the root universe
    if(present(n0)) then
      n = n0
      uni = self % universes(coords % lvl(n) % uniIdx) !*** BRAKES FORTRAN STANDARD. uni is undefined after execution
    else
      n = 1
      uni = self % rootUniverse          !*** BRAKES FORTRAN STANDARD. uni is undefined after execution
      call coords % takeAboveGeom()
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
        coords % regionID = c % uniqueID(instance)
        return

      ! Something has gone wrong...
      else
        print *, c % name()
        print *, c % latIdx()
        print *, c % uniIdx()
        print *, c % geometryIdx()
        print *, c % fillType()
        call fatalError('whichCell, geometry','Could not find the cell')
      end if
    end do

  end function whichCell

  !!
  !! Finds all the unique cell instances in the geometry and assigns them unique IDs
  !! This routine is used in initialisation
  !!
  subroutine ennumerateRegions(self)
    class(geometry), intent(inout)               :: self
    integer(shortInt)                            :: nRegions
    integer(shortInt), dimension(:), allocatable :: cellInstances
    integer(shortInt)                            :: cellIdx, latIdx, uniIdx, found
    type(coordList), dimension(:), allocatable   :: location
    type(coordList)                              :: searchLocation
    integer(shortInt)                            :: ID, i, n

    ! Track total number of regions
    nRegions = 0

    ! Instances of each cell (remains 0 for non-material cells)
    allocate(cellInstances(self % numCells))
    cellInstances(:) = 0

    ! Check each cell in the geometry by searching through all universe, starting with the root
    ! If filled by a material, enter the count subroutine
    ! Otherwise, continue
    print *,'Counting cell instances'
    do i=1,self % numCells
      uniIdx = self % rootUniverse % geometryIdx()
      latIdx = 0
      if (self % cells(i) % fillType == materialFill) then
        n = 0
        cellIdx = self % cells(i) % geometryIdx
        call self % countCell(cellIdx, uniIdx, latIdx, n)
        cellInstances(i) = n
        nRegions = nRegions + n
      end if
    end do

    self % numRegions = nRegions

    ! After counting each cell, must proceed to find them - provide with a coordList location
    print *,'Setting unique cell IDs and locations'
    call searchLocation % init([ZERO,ZERO,ZERO],[ONE,ZERO,ZERO])
    searchLocation % lvl(1) % uniIdx = self % rootUniverse % geometryIdx()
    ID = 1
    do i=1,self % numCells
      n = cellInstances(i)
      if (n > 1) then
        call searchLocation % takeAboveGeom()
        allocate(location(n))
        uniIdx = self % rootUniverse % geometryIdx()
        latIdx = 0
        cellIdx = self % cells(i) % geometryIdx
        found = 0
        call self % locateCell(cellIdx, uniIdx, latIdx, searchLocation, location, n, found)
        call self % cells(i) % setInstances(n,location,ID)
        ID = ID + n
        deallocate(location)
      end if
    end do

  end subroutine ennumerateRegions

  !!
  !! Count all the instances of a cell recursively
  !!
  recursive subroutine countCell(self, cellIdx, uniIdx, latIdx, n)
    class(geometry), intent(in)      :: self
    integer(shortInt), intent(in)    :: cellIdx
    integer(shortInt), intent(in)    :: uniIdx
    integer(shortInt), intent(in)    :: latIdx
    integer(shortInt), intent(inout) :: n
    integer(shortInt)                :: newUniIdx, newLatIdx, i, j, k, l
    type(cell_ptr)                   :: c

    do i=1,self % universes(uniIdx) % numCells
      c = self % universes(uniIdx) % cells(i)

      ! Material has been found - check against cellIdx under investigation
      if (c % fillType() == materialFill) then
        if (c % geometryIdx() == cellIdx) n = n + 1

      ! Descend to a lower universe
      else if (c % fillType() == universeFill) then
        newUniIdx = c % uniIdx()
        call self % countCell(cellIdx, newUniIdx, latIdx, n)

      ! Search through each universe of a lattice
      else if(c % fillType() == latticeFill) then
        newLatIdx = c % latIdx()
        do j = 1,self % lattices(newLatIdx) % extent(1)
          do k = 1,self % lattices(newLatIdx) % extent(2)
            do l = 1,self % lattices(newLatIdx) % extent(3)
              newUniIdx = self % lattices(newLatIdx) % universes(j,k,l) % geometryIdx()
              call self % countCell(cellIdx, newUniIdx, newLatIdx, n)
            end do
          end do
        end do

      ! An invalid fill has been found - something is terribly wrong
      else if (c % fillType() /= outsideFill) then
        call fatalError('countCell, geometry','Cannot identify cell contents')
      end if
    end do
    call c % kill()

  end subroutine countCell

  !!
  !! Locate all cell instances
  !!
  recursive subroutine locateCell(self, cellIdx, uniIdx, latIdx, searchLocation, location, n, found)
    class(geometry), intent(in)                  :: self
    integer(shortInt), intent(in)                :: cellIdx
    integer(shortInt), intent(in)                :: uniIdx
    integer(shortInt), intent(in)                :: latIdx
    type(coordList), intent(inout)               :: searchLocation
    type(coordList), dimension(:), intent(inout) :: location
    integer(shortInt), intent(in)                :: n
    integer(shortInt), intent(inout)             :: found
    integer(shortInt)                            :: i, j, k, l, nesting, &
                                                    ijkIdx, newLatIdx, newUniIdx
    type(cell_ptr)                               :: c

    do i=1,self % universes(uniIdx) % numCells
      c = self % universes(uniIdx) % cells(i)
      nesting = searchLocation % nesting
      searchLocation % lvl(nesting) % cellIdx = c % geometryIdx()

      ! Material has been found - check against cellIdx under investigation
      if (c % fillType() == materialFill) then
        if (c % geometryIdx() == cellIdx) then
          found = found + 1
          location(found) = searchLocation
        end if

      ! Descend to a lower universe
      else if (c % fillType() == universeFill) then
        newUniIdx = c % uniIdx()
        call searchLocation % addLevel([ZERO,ZERO,ZERO], newUniIdx)
        call self % locateCell(cellIdx, newUniIdx, latIdx, searchLocation, location, n, found)

      ! Search through each universe of a lattice
      else if(c % fillType() == latticeFill) then
        newLatIdx = c % latIdx()
        latSearch: do j = 1,self % lattices(newLatIdx) % extent(1)
          do k = 1,self % lattices(newLatIdx) % extent(2)
            do l = 1,self % lattices(newLatIdx) % extent(3)
              ijkIdx = self % lattices(newLatIdx) % getijkIdx([j,k,l])
              newUniIdx = self % lattices(newLatIdx) % universes(j,k,l) % geometryIdx()
              call searchLocation % addLevel([ZERO,ZERO,ZERO], newUniIdx, newlatIdx, ijkIdx)
              call self % locateCell(cellIdx, newUniIdx, newLatIdx, searchLocation, location, n, found)
              if (found == n) exit latSearch
            end do
          end do
        end do latSearch

      ! An invalid fill has been found - something is terribly wrong
      else if (c % fillType() /= outsideFill) then
        call fatalError('locateCell, geometry','Cannot identify cell contents')
      end if

      ! All instances have been found - exit the search!
      if (found == n) exit

    end do
    call c % kill()

    ! The current universe on this level has been fully explored - ascend a nesting level
    nesting = searchLocation % nesting
    if (nesting > 1) call searchLocation % decreaseNesting(nesting - 1)

  end subroutine locateCell

  !!
  !! Assumes that the geometry is bounded by a cuboid
  !! Details of boundary box provided by user
  !!
  subroutine constructBoundingBox(self, orig, a)
    class(geometry), intent(inout)          :: self
    real(defReal), intent(in), dimension(3) :: orig
    real(defReal), intent(in), dimension(3) :: a
    call self % boundingBox % init(orig, a)
    self % boxDefined = .TRUE.
  end subroutine constructBoundingBox

  !!
  !! Calculate cell volumes by sampling numPoints points in the geometry
  !! subject to the geometry having a bounding box
  !!
  !! Will eventually require parallelisation
  !!
  subroutine calculateVolumes(self, numPoints, seed)
    class(geometry), intent(inout)              :: self
    integer(shortInt), intent(in)               :: numPoints
    integer(longInt), intent(in), optional      :: seed
    integer(longInt)                            :: s
    real(defReal), dimension(3)                 :: minPoint, width, testPoint, dummyU
    integer(longInt), dimension(:), allocatable :: score
    real(defReal), dimension(:), allocatable    :: cellVols
    integer(shortInt)                           :: i, j, regionIdx
    type(rng)                                   :: rand
    type(box)                                   :: bb
    type(cell_ptr)                              :: c
    type(coordList)                             :: coords

    if(.NOT. self % boxDefined) call fatalError('calculateVolumes, geometry',&
                                'No bounding box has been defined')

    allocate(cellVols(self % numRegions))
    allocate(score(self % numRegions))
    score = 0
    cellVols = ZERO

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

      ! Return a base cell which either contains a material or sits outside the geometry
      c = self % whichCell(coords)

      if (c % insideGeom()) then
        regionIdx = coords % regionID
        score(regionIdx) = score(regionIdx) + 1
      end if
    end do
    call c % kill()

    cellVols = ONE* score * width(1) * width(2) * width(3) / numPoints

    ! Assign volumes to each material cell instance
    ! Loop through cells to find which have material fills
    ! If so, loop through each instance, checking their regionID
    ! and providing the corresponding volume
    do i = 1,self%numCells
      if (self % cells(i) % fillType == materialFill) then
        do j = 1,self % cells(i) % instances
          self % cells(i) % volume(j) = cellVols(self % cells(i) % uniqueID(j))
        end do
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
          colourMatrix(i,j) = coords % matIdx
        end if
      end do
    end do

    call c % kill()

  end function slicePlot

  !!
  !! Obtain a 3D array for use in voxel plotting
  !! Obtains material coloured voxel originating from a rear-bottom-left
  !! corner of the geometry and going forward by width in each direction
  !! Populates an outputMesh instance with material values
  !!
  function voxelPlot(self, nVox, corner, width) result(colourMatrix)
    class(geometry), intent(in)                      :: self
    integer(shortInt), dimension(3), intent(in)      :: nVox
    real(defReal), dimension(3), intent(in)          :: corner
    real(defReal), dimension(3), intent(in)          :: width
    integer(shortInt), dimension(:,:,:), allocatable :: colourMatrix
    real(defReal), dimension(3)                      :: x0, x, u
    integer(shortInt)                                :: i, j, k
    type(coordList)                                  :: coords
    type(cell_ptr)                                   :: c

    allocate(colourMatrix(nVox(1),nVox(2),nVox(3)))

    ! Create the co-ordinate points and widths
    x0 = corner + width*HALF
    x = x0
    u = [SQRT2_2, SQRT2_2, ZERO]

    call coords % init(x, u)

    ! Loop over each direction
    do k = 1,nVox(3)
      do j = 1,nVox(2)
        do i = 1,nVox(1)
          x = [x0(1) + (i-1)*width(1), x0(2) + (j-1)*width(2), x0(3) + (k-1)*width(3)]
          call coords % init(x, u)
          c = self % whichCell(coords)
          if (c % insideGeom()) then
            colourMatrix(i,j,k) = coords % matIdx
          else
            colourMatrix(i,j,k) = 0
          end if
        end do
      end do
    end do

    call c % kill()

  end function voxelPlot

  !!
  !! Given a dictionary describing a cell, construct a cell for the cell array
  !!
  subroutine cellFromDict(self, dict, geometryIdx, name)
    class(geometry), intent(inout)                :: self
    class(dictionary), intent(in)                 :: dict
    integer(shortInt), intent(in)                 :: geometryIdx
    character(*), intent(in)                      :: name
    integer(shortInt), dimension(:), allocatable  :: surfaces
    logical(defBool), dimension(:), allocatable   :: halfspaces
    class(surface_ptr), dimension(:), allocatable :: surfArray
    integer(shortInt)                             :: id, testID, i, j, fillType
    logical(defBool)                              :: insideGeom
    character(nameLen)                            :: fillName

    ! Determine how many surfaces the cell contains and set their halfspace
    allocate(surfaces(size(dict % getIntArray('surfaces'))))
    surfaces = dict % getIntArray('surfaces')
    allocate(surfArray(size(surfaces)))
    allocate(halfspaces(size(surfaces)))
    do i=1,size(surfaces)
      if(surfaces(i) > 0)then     ! *** Remove the loop. Use array expr. instead
        halfspaces(i) = .TRUE.
      else if(surfaces(i) < 0)then
        halfspaces(i) = .FALSE.
      else
        call fatalError('cellFromDict, geometry',&
        'Cannot have a surface with 0 halfspace in cell definition')
      end if
      do j=1,self % numSurfaces
        testID = self % surfaces(j) % id()
        if(testID == abs(surfaces(i))) surfArray(i) = self % surfaces(j)
      end do
    end do

    id = dict % getInt('id')
    fillName = dict % getChar('filltype')
    insideGeom = .TRUE.
    if(fillName == 'mat') then
      fillType = materialFill
    else if(fillName == 'uni') then
      fillType = universeFill
    else if(fillName == 'lat') then
      fillType = latticeFill
    else if(fillName == 'outside') then
      fillType = outsideFill
    else
      call fatalError('cellFromDict, geometry','Invalid fill type provided')
    end if

    ! Initialise the cell without filling: fill after all other geometry construction completed
    call self % cells(geometryIdx) % init(surfArray, halfspaces, id, fillType, geometryIdx, name)

    do i=1,size(surfaces)
      call surfArray(i) % kill()
    end do

  end subroutine cellFromDict

  !!
  !! Given a dictionary describing a universe, construct a universe for the universe array
  !!
  subroutine universeFromDict(self, dict, geometryIdx, name, isRoot)
    class(geometry), intent(inout)               :: self
    class(dictionary), intent(in)                :: dict
    integer(shortInt), intent(in)                :: geometryIdx
    character(*), intent(in)                     :: name
    logical(defBool), intent(in)                 :: isRoot
    integer(shortInt)                            :: id
    integer(shortInt)                            :: numCells, i, j, cellID, testId
    integer(shortInt), dimension(:), allocatable :: cellIDArray
    class(cell_ptr), dimension(:), allocatable   :: cellPtrArray
    logical(defBool)                             :: foundCell, foundOutside

    ! Determine how many cells the universe contains
    numCells = size(dict % getIntArray('cells'))
    allocate(cellIDArray(numCells))
    allocate(cellPtrArray(numCells))
    cellIDArray = dict % getIntArray('cells')

    ! Find all cells which occupy the universe
    do i = 1,numCells
      cellID = cellIDArray(i)
      foundCell = .FALSE.
      do j=1,self % numCells
        testId = self % cells(j) % id
        if(testId == cellID) then
          cellPtrArray(i) = self % cells(j)
          foundCell = .TRUE.
          exit
        end if
      end do
      if(.NOT.foundCell)then
        call fatalError('universeFromDict, geometry','Could not find a cell contained by a universe')
      end if
    end do

    ! If this universe is the root, make sure it contains at least two cells
    ! One of these cells must be the outside cell
    if (isRoot) then
      if (numCells < 2) call fatalError('universeFromDict, geometry',&
                                        'The base universe must contain at least two cells')
      foundOutside = .FALSE.
      do i =1,numCells
        if (cellPtrArray(i) % fillType() == outsideFill ) foundOutside = .TRUE.
      end do
      if (.NOT. foundOutside) call fatalError('universeFromDict, geometry',&
                                              'The base universe must contain the outside cell')
    end if

    ! Read unique universe ID
    id = dict % getInt('id')

    ! Does not require an input offset by default
    if(dict % isPresent('offset')) then
      call self % universes(geometryIdx) % &
      init(dict % getRealArray('offset'), cellPtrArray, id, geometryIdx, isRoot, name)
    else
      call self % universes(geometryIdx) % &
      init([ZERO,ZERO,ZERO], cellPtrArray, id, geometryIdx, isRoot, name)
    end if

    do i=1,numCells
      call cellPtrArray(i) % kill()
    end do

  end subroutine universeFromDict

  !!
  !! Given a dictionary describing a lattice, construct a lattice for the lattice array
  !!
  subroutine latticeFromDict(self, dict, name, geometryIdx)
    class(geometry), intent(inout)                     :: self
    class(dictionary), intent(in)                      :: dict
    character(*), intent(in)                           :: name
    integer(shortInt), intent(in)                      :: geometryIdx
    integer(shortInt), dimension(3)                    :: latShape
    integer(shortInt), dimension(:), allocatable       :: latIDArray
    class(universe_ptr), dimension(:,:,:), allocatable :: latUniverses
    integer(shortInt)                                  :: i,j,k,l,id, testID
    logical(defBool)                                   :: foundUni, is3D
    real(defReal), dimension(3)                        :: pitch, corner

    ! Determine the lattice structure
    latShape = dict % getIntArray('shape')
    if(any(latShape < 1)) call fatalError('latFromDict, geometry','Invalid lattice shape provided')
    allocate(latUniverses(latShape(1),latShape(2),latShape(3)))
    allocate(latIDArray(latShape(1)*latShape(2)*latShape(3)))
    latIDArray = dict % getIntArray('universes')

    if (latShape(3) > 1) then
      is3D = .TRUE.
    else
      is3D = .FALSE.
    end if
    ! Populate the universe pointer array
    do k=1,latShape(3)
      do i=1,latShape(1)
        do j=1,latShape(2)
          id = latIDArray(j + (i-1)*latShape(2) + (k-1)*latShape(1)*latShape(2))
          ! Search universe array to find ID corresponding to latUniverseID
          foundUni = .FALSE.
          do l=1,self % numUniverses
            testID = self % universes(l) % id
            if (testID == id) then
              latUniverses(i,j,k) = self % universes(l)
              foundUni = .TRUE.
              exit
            end if
          end do
          if (.NOT.foundUni) call &
          fatalError('latFromDict, geometry','Could not find a universe in the lattice')
        end do
      end do
    end do

    pitch = dict % getRealArray('pitch')
    corner = dict % getRealArray('corner')
    id = dict % getInt('id')
    call self % lattices(geometryIdx) % init(pitch, corner, latUniverses, id, is3D, name)
    do i=1,latShape(1)
      do j=1,latShape(2)
        do k=1,latShape(3)
          call latUniverses(i,j,k) % kill()
        end do
      end do
    end do

  end subroutine latticeFromDict

  !!
  !! Fills cells in the cell array with a universe
  !!
  subroutine fillCellUni(self, geometryIdx, dict)
    class(geometry), intent(inout)            :: self
    integer(shortInt), intent(in)             :: geometryIdx
    class(dictionary), intent(in)             :: dict
    integer(shortInt)                         :: i, id, testID

    id = dict % getInt('universe')
    do i = 1,self % numUniverses
      testID = self % universes(i) % id
      if (id == testID) then
        self % cells(geometryIdx) % uniIdx = i
        return
      end if
    end do

  end subroutine fillCellUni

  !!
  !! Fills cells in the cell array with a lattice
  !!
  subroutine fillCellLat(self, geometryIdx, dict)
    class(geometry), intent(inout)           :: self
    integer(shortInt), intent(in)            :: geometryIdx
    class(dictionary), intent(in)            :: dict
    integer(shortInt)                        :: i, id, testID

    id = dict % getInt('lattice')
    do i = 1,self % numLattices
      testID = self % lattices(i) % id
      if (id == testID) then
        self % cells(geometryIdx) % latIdx = i
        return
      end if
    end do

  end subroutine fillCellLat

  !!
  !! Applies the boundary conditions to the bounding surfaces as defined in the BC dictionary
  !!
  subroutine setBoundaryConditions(self, geometryIdx, dict)
    class(geometry), intent(inout)  :: self
    integer(shortInt), intent(in)   :: geometryIdx
    class(dictionary), intent(in)   :: dict
    integer(shortInt)               :: i,j
    integer(shortInt), dimension(6) :: BC

    ! Identify the boundary conditions
    BC= dict % getIntArray('boundary')
    if(self % cells(geometryIdx) % numSurfaces > 1) &
    call fatalError('applyBCs','Outside cell must only be composed of one surface')
    call self % cells(geometryIdx) % surfaces(1) % setBoundaryConditions(BC)
    self % boundarySurface = self % cells(geometryIdx) % surfaces(1)

  end subroutine setBoundaryConditions

  !!
  !! Fill cells containing a material
  !!
  subroutine fillMaterialCells(self, cellDict, materials)
    class(geometry), intent(inout)                :: self
    class(dictionary), intent(in)                 :: cellDict
    class(nuclearData), intent(inout)             :: materials
    character(nameLen), dimension(:), allocatable :: keys
    type(dictionary)                              :: subDict
    character(nameLen)                            :: name
    integer(shortInt)                             :: idx, i


    call cellDict % keysDict(keys)

    do i=1,self%numCells
      if (self % cells(i) % fillType == materialFill) then
        subDict = cellDict % getDict(keys(i))
        name = subDict % getChar('mat')
        idx = materials % getIdx(name)
        call self % cells(i) % fill(idx)
        !call materials % addActiveMaterial(name)
      end if
    end do
    deallocate(keys)
    call subDict % kill()

  end subroutine fillMaterialCells

end module geometry_class
