!
! Module for initialising arrays to be fed to geometry
!
module initialiseGeometryStructures
  use genericProcedures
  use numPrecision
  use universalVariables

  use surface_class
  use box_class
  use cylinder_class
  use sphere_class
  use plane_class
  use truncatedCylinder_class
  use squareCylinder_class

  use cell_class
  use universe_class
  use lattice_class

  use geometry_class

  use nuclearData_inter, only : nuclearData

  use dictionary_class

  implicit none
  private

  public :: initGeometryFromDict

contains

  !!
  !! Initialise the geometry on provision of a dictionary
  !!
  subroutine initGeometryFromDict(geom, surfaces, cells, universes, lattices, max_nest, rootIdx, dict, materials)
    class(geometry), intent(inout)                             :: geom
    class(dictionary), intent(inout)                           :: dict
    type(dictionary)                                           :: surfDict, cellDict, uniDict,&
                                                                  latDict, bcDict, uniDaughterDict
    character(100), dimension(:), allocatable                  :: keys
    integer(shortInt)                                          :: i,j, id, testID, &
                                                                  outsideDefined, numSurfaces, &
                                                                  numCells, numUniverses, &
                                                                  numLattices, fillType, outsideIdx
    logical(defBool)                                           :: foundRoot
    class(surface_ptr), dimension(:), allocatable, intent(out) :: surfaces
    class(cell), dimension(:), allocatable, intent(out)        :: cells
    class(universe), dimension(:), allocatable, intent(out)    :: universes
    class(lattice), dimension(:), allocatable, intent(out)     :: lattices
    class(nuclearData), intent(inout)                          :: materials
    integer(shortInt), intent(out)                             :: max_nest, rootIdx

    print *,'Reading geometry from dictionary'

    ! Read how many surfaces exist
    surfDict = dict % getDict('surfaces')
    numSurfaces = size(surfDict % keysDict())
    allocate(surfaces(numSurfaces))
    allocate(keys(numSurfaces))
    keys = surfDict % keysDict()

    ! Construct surfaces and store in the surface array
    print *,'Constructing ',numSurfaces,' surfaces'
    do i=1,numSurfaces
      call surfaceFromDict(surfDict % getDict(keys(i)), keys(i), surfaces(i))
    end do
    deallocate(keys)
    call surfDict % kill()

    ! Ensure no surfaces have the same id's
    if (numSurfaces>1) then
      do i=1,numSurfaces-1
        id = surfaces(i) % id()
        do j=i+1,numSurfaces
          testId = surfaces(j) % id()
          if(testId==id) call fatalError('initGeomArraysDict','Surfaces defined with identical ids')
        end do
      end do
    end if

    ! Read how many cells exist
    cellDict = dict % getDict('cells')
    numCells = size(cellDict % keysDict())
    allocate(cells(numCells))
    allocate(keys(numCells))
    keys = cellDict % keysDict()

    ! Construct cells and store in the cell array
    ! Also ensure that the outside cell is defined
    outsideDefined = 0
    print *,'Constructing cells'
    do i=1,numCells
      call cellFromDict(surfaces, cellDict % getDict(keys(i)), i, keys(i), cells(i))
      if (cells(i) % fillType == outsideFill) then
        outsideDefined = outsideDefined + 1
      end if
    end do
    deallocate(keys)

    if (outsideDefined < 1) then
      call fatalError('initGeomArraysDict','Must define one outside cell')
    else if (outsideDefined > 1) then
      call fatalError('initGeomArraysDict','More than one outside cell defined')
    end if

    if (numCells < 2) then
      call fatalError('initGeomArraysDict','Must define at least an outside and a material cell')
    end if

    ! Delay termination of cellDict termination until cells are filled with materials/universes/lattices
    ! Ensure no cells have the same id's
    do i=1,numCells-1
      id = cells(i) % id
      do j=i+1,numCells
        testId = cells(j) % id
        if(testId==id) call fatalError('initGeomArraysDict','Cells defined with identical ids')
      end do
    end do

    ! Read how many universes exist
    uniDict = dict % getDict('universes')
    numUniverses = size(uniDict % keysDict())
    allocate(universes(numUniverses))
    allocate(keys(numUniverses))
    keys = uniDict % keysDict()

    ! Construct universes and store in the universe array
    print *,'Constructing universes'
    foundRoot = .FALSE.
    do i=1,numUniverses
      uniDaughterDict = uniDict % getDict(keys(i))
      if(uniDaughterDict % isPresent('root')) then
        if (uniDaughterDict % getInt('root') == 1) then
          call universeFromDict(cells, uniDaughterDict, i, keys(i), .TRUE., universes(i))
          rootIdx = i
          if (.not.foundRoot) then
            foundRoot = .TRUE.
          else
            call fatalError('initGeomArraysDict','More than one root universe defined')
          end if
        end if
      else
        call universeFromDict(cells, uniDaughterDict, i, keys(i), .FALSE., universes(i))
      end if
    end do
    deallocate(keys)
    call uniDict % kill()

    ! Ensure that a root universe has been defined
    if(.not.foundRoot)then
      call fatalError('initGeomArraysDict','The root universe has not been defined')
    end if

    ! Ensure no universes have the same id's
    if (numUniverses>1) then
      do i=1,numUniverses-1
        id = universes(i) % id
        do j=i+1,numUniverses
          testId = universes(j) % id
          if(testId==id) call fatalError('initGeomArraysDict','Universes defined with identical ids')
        end do
      end do
    end if

    ! Check if lattices are included in the geometry
    if(dict % isPresent('lattices')) then
      print *,'Constructing lattices'

      ! Read how many lattices exist
      latDict = dict % getDict('lattices')
      numLattices = size(latDict % keysDict())
      allocate(lattices(numLattices))
      allocate(keys(numLattices))
      keys = latDict % keysDict()

      ! Construct lattices and store in the lattice array
      do i=1,numLattices
        call latticeFromDict(universes, latDict % getDict(keys(i)), keys(i), lattices(i))
      end do
      deallocate(keys)
      call latDict % kill()

      ! Ensure no universes have the same id's
      if (numLattices > 1) then
        do i=1,numLattices-1
          id = lattices(i) % id
          do j=i+1,numLattices
            testId = lattices(j) % id
            if(testId==id) call fatalError('initGeomArraysDict','Lattices defined with identical ids')
          end do
        end do
      end if
    else
      print *,'No lattices in the geometry'
    end if

    ! Fill all cells given that universes, lattices, and materials are all defined
    allocate(keys(numCells))
    keys = cellDict % keysDict()
    print *,'Filling cells with universes and lattices'
    do i=1,numCells
      fillType = cells(i) % fillType
      if (fillType == universeFill) then
        call fillCellUni(cells(i), cellDict % getDict(keys(i)), universes)
      else if (fillType == latticeFill) then
        call fillCellLat(cells(i), cellDict % getDict(keys(i)), lattices)
      else if (fillType == materialFill) then
        cycle ! Fill with materials after number of cell instances have been identified
        !call fillCellMat(cells(i), cellDict % getDict(keys(i)), materials)
      else if (fillType == outsideFill) then
        ! Read and apply boundary conditions
        print *,'Applying boundary conditions'
        call setBoundaryConditions(cells(i), dict)
        outsideIdx = i
        cycle
      else
        call fatalError('fillCell, initialiseGeometryStructures','Cell has an incorrect fill type')
      end if
    end do

    ! Calculate maximum nesting
    print *,'Calculating maximum geometry nesting'
    max_nest = hardcoded_max_nest

    ! *** CalcMaxNesting is pending debbuging. May appear at some point
    !call calcMaxNesting(max_nest,universes,lattices, rootInd)
    print *,'Geometry contains ',max_nest,' nesting levels'

    ! Initialise gometry
    call geom % init(surfaces, cells, universes, rootIdx, cells(outsideIdx) % surfaces(1), lattices)

    ! Identify all unique base cell instances
    ! Required for regional tallies, MOC flux regions, and eventual burn-up

    ! Fill cells with their materials
    call cellDict % kill()
    call dict % kill()
    deallocate(keys)

  end subroutine initGeometryFromDict

  !!
  !! Given a dictionary describing a surface, construct a surface for the surrface array
  !!
  subroutine surfaceFromDict(dict, name, surf)
    class(dictionary), intent(in)     :: dict
    character(*), intent(in)          :: name
    integer(shortInt)                 :: id
    type(box), pointer                :: boxObj
    type(sphere), pointer             :: sphereObj
    type(xCylinder), pointer          :: xCylObj
    type(yCylinder), pointer          :: yCylObj
    type(zCylinder), pointer          :: zCylObj
    type(xSquareCylinder), pointer    :: xSquCylObj
    type(ySquareCylinder), pointer    :: ySquCylObj
    type(zSquareCylinder), pointer    :: zSquCylObj
    type(xTruncatedCylinder), pointer :: xTruncCylObj
    type(yTruncatedCylinder), pointer :: yTruncCylObj
    type(zTruncatedCylinder), pointer :: zTruncCylObj
    type(xPlane), pointer             :: xPlaneObj
    type(yPlane), pointer             :: yPlaneObj
    type(zPlane), pointer             :: zPlaneObj
    type(plane), pointer              :: planeObj
    character(100)                    :: surfType
    real(defReal), dimension(4)       :: coeff
    real(defReal), dimension(3)       :: origin, halfwidth
    real(defReal)                     :: radius, halfheight, x, y, z
    type(surface_ptr), intent(inout)  :: surf

    surfType = dict % getChar('type')
    id = dict % getInt('id')
    if(id<1) call fatalError('surfaceFromDict, initialiseGeometryStructures','Invalid surface id provided')


    if (surfType == 'box') then
      allocate(boxObj)
      halfwidth = dict % getRealArray('halfwidth')
      origin = dict % getRealArray('origin')
      call boxObj % init(origin, halfwidth, id, name)
      surf = boxObj
      return

    else if (surfType == 'xSquareCylinder') then
      allocate(xSquCylObj)
      halfwidth = dict % getRealArray('halfwidth')
      origin = dict % getRealArray('origin')
      call xSquCylObj % init(origin, halfwidth, id, name)
      surf = xSquCylObj
      return

    else if (surfType == 'ySquareCylinder') then
      allocate(ySquCylObj)
      halfwidth = dict % getRealArray('halfwidth')
      origin = dict % getRealArray('origin')
      call ySquCylObj % init(origin, halfwidth, id, name)
      surf = ySquCylObj
      return

    else if (surfType == 'zSquareCylinder') then
      allocate(zSquCylObj)
      halfwidth = dict % getRealArray('halfwidth')
      origin = dict % getRealArray('origin')
      call zSquCylObj % init(origin, halfwidth, id, name)
      surf = zSquCylObj
      return

    else if (surfType == 'sphere') then
      allocate(sphereObj)
      radius = dict % getReal('radius')
      origin = dict % getRealArray('origin')
      call sphereObj % init(origin, radius, id, name)
      surf = sphereObj
      return

    else if (surfType == 'xCylinder') then
      allocate(xCylObj)
      radius = dict % getReal('radius')
      origin = dict % getRealArray('origin')
      call xCylObj % init(radius, origin, id, name)
      surf = xCylObj
      return

    else if (surfType == 'yCylinder') then
      allocate(yCylObj)
      radius = dict % getReal('radius')
      origin = dict % getRealArray('origin')
      call yCylObj % init(radius, origin, id, name)
      surf = yCylObj
      return

    else if (surfType == 'zCylinder') then
      allocate(zCylObj)
      radius = dict % getReal('radius')
      origin = dict % getRealArray('origin')
      call zCylObj % init(radius, origin, id, name)
      surf = zCylObj
      return

    else if (surfType == 'xTruncCylinder') then
      allocate(xTruncCylObj)
      radius = dict % getReal('radius')
      halfheight = dict % getReal('halfheight')
      origin = dict % getRealArray('origin')
      call xTruncCylObj % init(origin, halfheight, radius, id, name)
      surf = xTruncCylObj
      return

    else if (surfType == 'yTruncCylinder') then
      allocate(yTruncCylObj)
      radius = dict % getReal('radius')
      halfheight = dict % getReal('halfheight')
      origin = dict % getRealArray('origin')
      call yTruncCylObj % init(origin, halfheight, radius, id, name)
      surf = yTruncCylObj
      return

    else if (surfType == 'zTruncCylinder') then
      allocate(zTruncCylObj)
      radius = dict % getReal('radius')
      halfheight = dict % getReal('halfheight')
      origin = dict % getRealArray('origin')
      call zTruncCylObj % init(origin, halfheight, radius, id, name)
      surf = zTruncCylObj
      return

    else if (surfType == 'xPlane') then
      allocate(xPlaneObj)
      x = dict % getReal('x')
      call xPlaneObj % init(x, id, name)
      surf = xPlaneObj
      return

    else if (surfType == 'yPlane') then
      allocate(yPlaneObj)
      y = dict % getReal('y')
      call yPlaneObj % init(y, id, name)
      surf = yPlaneObj
      return

    else if (surfType == 'zPlane') then
      allocate(zPlaneObj)
      z = dict % getReal('z')
      call xPlaneObj % init(z, id, name)
      surf = zPlaneObj
      return

    else if (surfType == 'plane') then
      allocate(planeObj)
      coeff = dict % getRealArray('coeff')
      call planeObj % init(coeff, id, name)
      surf = planeObj
      return

    else
      call fatalError('surfaceFromDict','Did not recognise surface type')
    end if
  end subroutine surfaceFromDict

  !!
  !! Given a dictionary describing a cell, construct a cell for the cell array
  !!
  subroutine cellFromDict(surfacesGeom, dict, geometryIdx, name, c)
    class(dictionary), intent(in)                 :: dict
    integer(shortInt), intent(in)                 :: geometryIdx
    class(surface_ptr), dimension(:), intent(in)  :: surfacesGeom
    character(*), intent(in)                      :: name
    type(cell), intent(inout)                     :: c
    integer(shortInt), dimension(:), allocatable  :: surfaces
    logical(defBool), dimension(:), allocatable   :: halfspaces
    class(surface_ptr), dimension(:), allocatable :: surfArray
    integer(shortInt)                             :: id, testID, i, j, fillType
    logical(defBool)                              :: insideGeom
    character(100)                                :: fillName

    ! Determine how many surfaces the cell contains and set their halfspace
    allocate(surfaces(size(dict % getIntArray('surfaces'))))
    surfaces = dict % getIntArray('surfaces')
    allocate(surfArray(size(surfaces)))
    allocate(halfspaces(size(surfaces)))
    do i=1,size(surfaces)
      if(surfaces(i)>0)then
        halfspaces(i) = .TRUE.
      else if(surfaces(i)<0)then
        halfspaces(i) = .FALSE.
      else
        call fatalError('cellFromDict, initialiseGeometryStructures',&
        'Cannot have a surface with 0 halfspace in cell definition')
      end if
      do j=1,size(surfacesGeom)
        testID = surfacesGeom(j) % id()
        if(testID == abs(surfaces(i))) surfArray(i) = surfacesGeom(j)
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
      call fatalError('cellFromDict, initialiseGeometryStructures','Invalid fill type provided')
    end if

    ! Initialise the cell without filling: fill after all other geometry construction completed
    call c % init(surfArray, halfspaces, id, fillType, geometryIdx, name)

    do i=1,size(surfaces)
      call surfArray(i) % kill()
    end do

  end subroutine cellFromDict

  !!
  !! Given a dictionary describing a universe, construct a universe for the universe array
  !!
  subroutine universeFromDict(cellArray, dict, geometryIdx, name, isRoot, uni)
    class(cell), dimension(:)                    :: cellArray
    class(dictionary), intent(in)                :: dict
    integer(shortInt), intent(in)                :: geometryIdx
    character(*), intent(in)                     :: name
    logical(defBool), intent(in)                 :: isRoot
    integer(shortInt)                            :: id
    type(universe), intent(inout)                :: uni
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
    do i=1,numCells
      cellID = cellIDArray(i)
      foundCell = .FALSE.
      do j=1,size(cellArray)
        testId = cellArray(j) % id
        if(testId == cellID) then
          cellPtrArray(i) = cellArray(j)
          foundCell = .TRUE.
          exit
        end if
      end do
      if(.not.foundCell)then
        call fatalError('uniFromDict, geometry','Could not find a cell contained by a universe')
      end if
    end do

    ! If this universe is the root, make sure it contains at least two cells
    ! One of these cells must be the outside cell
    if (isRoot) then
      if (numCells < 2) call fatalError('uniFromDict','The base universe must contain at least two cells')
      foundOutside = .FALSE.
      do i =1,numCells
        if (cellPtrArray(i) % fillType() == outsideFill ) foundOutside = .TRUE.
      end do
      if (.NOT. foundOutside) call fatalError('uniFromDict','The base universe must contain the outside cell')
    end if

    ! Read unique universe ID
    id = dict % getInt('id')

    ! Does not require an input offset by default
    if(dict % isPresent('offset')) then
      call uni % init(dict % getRealArray('offset'), cellPtrArray, id, geometryIdx, isRoot, name)
    else
      call uni % init([ZERO,ZERO,ZERO], cellPtrArray, id, geometryIdx, isRoot, name)
    end if

    do i=1,numCells
      call cellPtrArray(i) % kill()
    end do

  end subroutine universeFromDict

  !!
  !! Given a dictionary describing a lattice, construct a lattice for the lattice array
  !!
  subroutine latticeFromDict(universeArray, dict, name, lat)
    class(universe), dimension(:), intent(in)          :: universeArray
    class(dictionary), intent(in)                      :: dict
    character(*), intent(in)                           :: name
    type(lattice), intent(inout)                       :: lat
    integer(shortInt), dimension(3)                    :: latShape
    integer(shortInt), dimension(:), allocatable       :: latIDArray
    class(universe_ptr), dimension(:,:,:), allocatable :: latUniverses
    integer(shortInt)                                  :: i,j,k,l,id, testID
    logical(defBool)                                   :: foundUni, is3D
    real(defReal), dimension(3)                        :: pitch, corner

    ! Determine the lattice structure
    latShape = dict % getIntArray('shape')
    if(any(latShape<1)) call fatalError('latFromDict, initialiseGeometryStructures','Invalid lattice shape provided')
    allocate(latUniverses(latShape(1),latShape(2),latShape(3)))
    allocate(latIDArray(latShape(1)*latShape(2)*latShape(3)))
    latIDArray = dict % getIntArray('universes')

    if (latShape(3)>1) then
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
          do l=1,size(universeArray)
            testID = universeArray(l) % id
            if (testID == id) then
              latUniverses(i,j,k) = universeArray(l)
              foundUni = .TRUE.
              exit
            end if
          end do
          if (.NOT.foundUni) call &
          fatalError('latFromDict, initialiseGeometryStructures',&
          'Could not find a universe in the lattice')
        end do
      end do
    end do

    pitch = dict % getRealArray('pitch')
    corner = dict % getRealArray('corner')
    id = dict % getInt('id')
    call lat % init(pitch, corner, latUniverses, id, is3D, name)
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
  subroutine fillCellUni(c, dict, universeArray)
    class(cell), intent(inout)                :: c
    class(dictionary), intent(in)             :: dict
    class(universe), dimension(:), intent(in) :: universeArray
    integer(shortInt)                         :: i, id, testID

    id = dict % getInt('universe')
    do i = 1,size(universeArray)
      testID = universeArray(i) % id
      if (id == testID) then
        c % uniIdx = i
        return
      end if
    end do

  end subroutine fillCellUni

  !!
  !! Fills cells in the cell array with a lattice
  !!
  subroutine fillCellLat(c, dict, latticeArray)
    class(cell), intent(inout)               :: c
    class(dictionary), intent(in)            :: dict
    class(lattice), dimension(:), intent(in) :: latticeArray
    integer(shortInt)                        :: i, id, testID

    id = dict % getInt('lattice')
    do i = 1,size(latticeArray)
      testID = latticeArray(i) % id
      if (id == testID) then
        c % latIdx = i
        return
      end if
    end do

  end subroutine fillCellLat

  !!
  !! Fills cells in the cell array with a material
  !!
  subroutine fillCellMat(c, dict,  materials)
    class(cell), intent(inout)        :: c
    class(dictionary), intent(in)     :: dict
    class(nuclearData), intent(inout) :: materials
    character(100)                    :: name
    integer(shortInt)                 :: idx

    name = dict % getChar('mat')
    idx = materials % getIdx(name)
    c % materialIdx = idx
    !call materials % addActiveMaterial(name)

  end subroutine fillCellMat

  !!
  !! Applies the boundary conditions to the bounding surfaces as defined in the BC dictionary
  !!
  subroutine setBoundaryConditions(outsideCell, dict)
    class(cell), intent(inout)      :: outsideCell
    class(dictionary), intent(in)   :: dict
    integer(shortInt)               :: i,j
    integer(shortInt), dimension(6) :: BC

    ! Identify the boundary conditions
    BC= dict % getIntArray('boundary')
    if(outsideCell % numSurfaces > 1) &
    call fatalError('applyBCs','Outside cell must only be composed of one surface')
    call outsideCell % surfaces(1) % setBoundaryConditions(BC)

  end subroutine setBoundaryConditions

  !!
  !! Find the nesting of the geometry for using co-ord lists
  !! **** PENDING DEBBUGING, SEGMENTATION AT THE MOMENT ***
  !!
!  recursive subroutine calcMaxNesting(max_nest, universes, lattices, uniInd)
!    class(universe), dimension(:), intent(in) :: universes
!    class(lattice), dimension(:), intent(in) :: lattices
!    integer(shortInt), intent(inout) :: max_nest
!    integer(shortInt), intent(in) :: uniInd
!    type(universe_ptr) :: uni
!    type(cell_ptr) :: c
!    integer(shortInt), dimension(3) :: ijk
!    integer(shortInt) :: cInd, newInd, latInd, test_nest, max_nest0, fillType, i, j, k
!
!    ! Include a simple protection against infinite nesting
!    if (max_nest > 10) call fatalError('calcMaxNesting','Geometry is more deeply nested than anticipated!')
!    print *,'Error here', max_nest
!
!    uni = universes(uniInd)
!    max_nest = max_nest + 1
!    max_nest0 = max_nest
!    test_nest = max_nest
!    print *,'uni numcells =',uni%numCells()
!    do cInd = 1,uni%numCells()
!      print *, cInd
!      c = uni % cells(cInd)
!      print *, associated(c % ptr), c % fillType()
!      print *,c % name()
!      if (c % fillType() == universeFill) then
!        newInd = c % uniInd()
!        call calcMaxNesting(test_nest, universes, lattices, newInd)
!        max_nest = max(max_nest, test_nest)
!        test_nest = max_nest0
!      else if (c % fillType() == latticeFill) then
!        latInd = c % latInd()
!        ijk = lattices(latInd) % extent
!        do i=1,ijk(1)
!          do j=1,ijk(2)
!            do k=1,ijk(3)
!              newInd = lattices(latInd) % universes(i,j,k) % geometryInd()
!              call calcMaxNesting(test_nest, universes, lattices, newInd)
!              max_nest = max(max_nest, test_nest)
!              test_nest = max_nest0
!            end do
!          end do
!        end do
!      else if (c % fillType() == materialFill) then
!        print *,'next'
!        cycle
!      else if (c % fillType() == outsideFill) then
!        cycle
!      else
!        call fatalError('calcMaxNesting, intialiseGeometryStructures',&
!                        'Could not identify cell fill type')
!      end if
!    end do
!    call uni % kill()
!    call c % kill()
!
!  end subroutine calcMaxNesting


end module initialiseGeometryStructures
