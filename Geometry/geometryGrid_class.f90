module geometryGrid_class

  use numPrecision
  use universalVariables
  use genericProcedures,  only : fatalError, numToChar
  use coord_class,        only : coordList, coord
  use dictionary_class,   only : dictionary
  use charMap_class,      only : charMap
  use geometry_inter,     only : geometry, distCache
  use csg_class,          only : csg
  use universe_inter,     only : universe
  use latUniverse_class,  only : latUniverse
  use surface_inter,      only : surface

  ! Nuclear Data
  use materialMenu_mod,   only : mm_matTemp => matTemp ,&
                                 mm_matFile => matFile ,&
                                 mm_init    => init    ,&
                                 mm_kill    => kill,&
                                 mm_matName => matName

  implicit none
  private

  !!
  !! A very simplified geometry model consisting of a uniform grid. This simplicity means we don't
  !! need to worry about nesting levels or complex shapes making things in general much faster.
  !! Useful for many IMC benchmarks when we need to split each material region into discrete zones
  !! to accurately model a time-evolving temperature field, especially when using a large number of
  !! zones.
  !!
  !! Makes use of csg_class and copies of geometryStd_class functions for initialisation only, such
  !! that input files can be virtually identical to geometryStd.
  !!
  !! Sample Dictionary Input:
  !!   geometry {
  !!     type geometryGrid;
  !!     dimensions (10 10 10);
  !!     <csg_class definition>
  !!    }
  !!
  !! This sample input will build the CSG geometry as in geometryStd_class, then construct a simple
  !! grid geometry automatically by with 10x10x10 cells, with the material in each cell equal to the
  !! material at the central point of that grid cell in the original CSG geometry. Each grid cell
  !! has a new instance of each material, even if multiple grid cells are contained in a single CSG
  !! material. MatIdxs will line up with localID if there are no VOID regions, if void is present
  !! then this will become out of sync - however matName e.g. 'mat63' correctly corresponds to
  !! position 63 (can be converted to ijk position similar to get_ijk in latUniverse_class).
  !!
  !! Interface:
  !!   Geometry Interface
  !!
  type, public, extends(geometry) :: geometryGrid
    type(csg) :: geom
    integer(shortInt), dimension(:), allocatable     :: latSizeN
    real(defReal), dimension(3)                      :: latPitch
    real(defReal), dimension(3)                      :: corner
    integer(shortInt), dimension(:,:,:), allocatable :: mats
    real(defReal), dimension(6)                      :: geomBounds
    integer(shortInt), dimension(:), allocatable     :: boundary

  contains
    ! Superclass procedures
    procedure :: init
    procedure :: kill
    procedure :: placeCoord
    procedure :: whatIsAt
    procedure :: bounds
    procedure :: move_noCache
    procedure :: move_withCache
    procedure :: moveGlobal
    procedure :: teleport
    procedure :: activeMats

    ! Public procedure unique to this class
    procedure :: matBounds

    ! Private procedures
    procedure, private :: initGridOnly
    procedure, private :: explicitBC
    procedure, private :: csg_diveToMat
    procedure, private :: csg_placeCoord
    procedure, private :: csg_whatIsAt
  end type geometryGrid

contains

  !!
  !! Initialise geometry
  !!
  !! See geometry_inter for details
  !!
  subroutine init(self, dict, mats, silent)
    class(geometryGrid), intent(inout)     :: self
    class(dictionary), intent(in)          :: dict
    type(charMap), intent(in)              :: mats
    logical(defBool), optional, intent(in) :: silent
    logical(defBool)                       :: loud
    real(defReal), dimension(6)            :: bounds
    class(surface), pointer                :: surf
    real(defReal), dimension(3)            :: r
    type(dictionary)                       :: matDict, tempDict
    real(defReal)                          :: volume
    character(nameLen)                     :: gridOnly
    integer(shortInt)                      :: i, j, k, z, N, matIdx, uniqueID, idxCounter, voidCounter
    character(100), parameter :: Here = 'init (geometryGrid_class.f90)'

    ! Choose whether to display messages
    if (present(silent)) then
      loud = .not.silent
    else
      loud = .true.
    end if

    ! Get geometry discretisation
    call dict % get(self % latSizeN, 'dimensions')
    if (size(self % latSizeN) /= 3) call fatalError(Here, 'Dimensions must be of size 3')

    ! Allocate space for material indexes
    allocate(self % mats(self % latSizeN(1),self % latSizeN(2),self % latSizeN(3)))

    ! Get boundary conditions
    call dict % get(self % boundary, 'boundary')
    if (size(self % boundary) /= 6) call fatalError(Here, 'boundary should be an array of size 6')

    ! Determine whether to completely rebuild geometry or to just build a grid
    call dict % getOrDefault(gridOnly, 'gridOnly','n')
    if (gridOnly == 'y') then
      ! Switch to dedicated subroutine
      call self % initGridOnly(dict)
      return
    else if (gridOnly /= 'n') then
      call fatalError(Here, 'Unrecognised value for gridOnly setting. Should be y or n.')
    end if

    ! Build the representation using CSG geometry
    call self % geom % init(dict, mats, silent)

    if (loud) then
      print *, "/\/\ CONVERTING CSG GEOMETRY TO GRID GEOMETRY /\/\"
    end if

    ! Get geometry bounds
    surf => self % geom % surfs % getPtr(self % geom % borderIdx)
    bounds = surf % boundingBox()
    self % geomBounds = bounds

    do i = 1, 3
      self % latPitch(i) = (bounds(i+3) - bounds(i)) / self % latSizeN(i)
    end do
    self % corner = bounds(1:3)
    volume = product(self % latPitch)

    ! Initialise dictionary of materials for materialMenu_mod initialisation
    N = product(self % latSizeN)
    call matDict % init(1)

    ! Loop through each grid cell
    idxCounter  = 0
    voidCounter = 0
    do k = 1, self % latSizeN(3)

      ! Flip in z axis for consistency with latUniverse_class
      z = self % latSizeN(3) - k + 1

      do j = 1, self % latSizeN(2)
        do i = 1, self % latSizeN(1)

          ! Get material at cell centre
          r = self % corner + [i-HALF,j-HALF,z-HALF]*self % latPitch
          call self % csg_whatIsAt(matIdx, uniqueID, r)

          ! Don't create a new material for void regions
          if (matIdx == VOID_MAT) then
            self % mats(i,j,z) = VOID_MAT
            ! Count number of void regions so that matName will match up with lattice position
            ! => easier data processing after obtaining results
            voidCounter = voidCounter + 1
            cycle
          end if

          ! Next matIdx
          idxCounter = idxCounter + 1

          ! Store in dictionary of new materials
          call tempDict % init(3)
          call tempDict % store('temp', mm_matTemp(matIdx))
          call tempDict % store('volume', volume)
          call tempDict % store('xsFile', mm_matFile(matIdx))
          call matDict % store('mat'//numToChar(idxCounter+voidCounter), tempDict)
          call tempDict % kill()

          ! Store matIdx in material array
          self % mats(i,j,z) = idxCounter

        end do
      end do
    end do

    ! Kill current geometry and materials
    call self % geom % kill()
    call mm_kill()

    ! Initialise new materials
    call mm_init(matDict)

    ! Kill material dictionary
    call matDict % kill()

    ! Print finish line
    if (loud) then
      print *, "\/\/ FINISHED BUILDING GRID GEOMETRY \/\/"
      print *, repeat('<>', MAX_COL/2)
    end if

  end subroutine init

  !!
  !! Generate grid geometry without giving any consideration to materials as is done in the full
  !! subroutine above. Assigns a unique value to each grid cell so that self % whatIsAt(matIdx)
  !! will return useful value.
  !!
  !! Unlike in full subroutine above, "bounds" is required in input dictionary.
  !!
  subroutine initGridOnly(self, dict)
    class(geometryGrid), intent(inout)       :: self
    class(dictionary), intent(in)            :: dict
    real(defReal), dimension(:), allocatable :: bounds
    integer(shortInt)                        :: i, j, k, idxCounter

    ! Get grid bounds and calculate basic properties
    call dict % get(bounds, 'bounds')
    self % geomBounds = bounds
    do i = 1, 3
      self % latPitch(i) = (bounds(i+3) - bounds(i)) / self % latSizeN(i)
    end do
    self % corner = bounds(1:3)

    ! Assign a unique value to each grid cell
    idxCounter = 1
    do k = 1, self % latSizeN(3)
      do j = 1, self % latSizeN(2)
        do i = 1, self % latSizeN(1)
          self % mats(i,j,k) = idxCounter
          idxCounter = idxCounter + 1
        end do
      end do
    end do

  end subroutine initGridOnly

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(geometryGrid), intent(inout) :: self

    call self % geom % kill()
    deallocate(self % mats)

  end subroutine kill

  !!
  !! Place coordinate list into geometry
  !!
  !! See geometry_inter for details
  !!
  subroutine placeCoord(self, coords)
    class(geometryGrid), intent(in) :: self
    type(coordList), intent(inout) :: coords
    integer(shortInt)              :: matIdx, uniqueID
    character(100), parameter :: Here = 'placeCoord (geometryGrid_class.f90)'

    call self % whatIsAt(matIdx, uniqueID, coords % lvl(1) % r)
    coords % matIdx = matIdx

    ! Extra unnecessary info for coords % isPlaced to return true
    coords % uniqueID = matIdx
    coords % nesting  = 1

  end subroutine placeCoord

  !!
  !! Find material and unique cell at a given location
  !! Optional direction input is not used to nudge particles across cells. All moves in this class
  !! are increased by NUDGE so particles should never be exactly on a surface.
  !!
  !! See geometry_inter for details
  !!
  subroutine whatIsAt(self, matIdx, uniqueID, r, u)
    class(geometryGrid), intent(in)                   :: self
    integer(shortInt), intent(out)                    :: matIdx
    integer(shortInt), intent(out)                    :: uniqueID
    real(defReal), dimension(3), intent(in)           :: r
    real(defReal), dimension(3), optional, intent(in) :: u
    integer(shortInt), dimension(3)                   :: ijk

    ! Avoid rounding error for leaked particles very close to boundary
    if (any(r < self % corner) .or. any(r > self % geomBounds(4:6))) then
      matIdx = OUTSIDE_MAT
      uniqueID = matIdx
      return
    end if

    ! Determine ijk location of cell
    ijk = floor((r - self % corner) / self % latPitch) + 1

    ! Get matIdx from array
    matIdx = self % mats(ijk(1),ijk(2),ijk(3))

    ! UniqueID not needed, set equal to matIdx to avoid compiler warning
    uniqueID = matIdx

  end subroutine whatIsAt

  !!
  !! Return Axis Aligned Bounding Box encompassing the geometry
  !!
  !! See geometry_inter for details
  !!
  function bounds(self)
    class(geometryGrid), intent(in) :: self
    real(defReal), dimension(6)     :: bounds

    bounds = self % geomBounds

  end function bounds

  !!
  !! Given coordinates placed in the geometry move point through the geometry
  !!
  !! See geometry_inter for details
  !!
  !! Uses explicit BC
  !!
  subroutine move_noCache(self, coords, maxDist, event)
    class(geometryGrid), intent(in) :: self
    type(coordList), intent(inout) :: coords
    real(defReal), intent(inout)   :: maxDist
    integer(shortInt), intent(out) :: event
    real(defReal)                  :: dist
    real(defReal), dimension(3)    :: r, u, r_bar
    character(100), parameter :: Here = 'move (geometryGrid_class.f90)'

    ! Calculate distance to next cell crossing
    r = coords % lvl(1) % r
    u = coords % lvl(1) % dir
    r_bar = (r - self % corner) / self % latPitch ! Normalise position within grid
    r_bar = r_bar - floor(r_bar)                  ! Normalise position within cell
    r_bar = (HALF - r_bar + sign(HALF, u)) * self % latPitch
    dist = minval(r_bar / u)                      ! Which direction will result in crossing

    if (maxDist < dist) then ! Moves within cell
      call coords % moveGlobal(maxDist)
      event = COLL_EV
      maxDist = maxDist ! Left for explicitness and compiler

      ! Place coords back into geometry
      call self % placeCoord(coords)

    else ! Move to next cell, increased by NUDGE to avoid numerical issues
      call coords % moveGlobal(dist + NUDGE)
      maxDist = dist + NUDGE

      ! Set matIdx
      call self % placeCoord(coords)

      ! Apply boundary conditions if leaving geometry
      if (coords % matIdx == OUTSIDE_MAT) then 
        event = BOUNDARY_EV
        call self % explicitBC(coords)

      else
        ! Cell crossing within geometry - no BCs needed
        event = CROSS_EV

      end if

    end if

  end subroutine move_noCache

  !!
  !! Given coordinates placed in the geometry move point through the geometry
  !!
  !! See geometry_inter for details
  !!
  !! Uses explicit BC
  !!
  subroutine move_withCache(self, coords, maxDist, event, cache)
    class(geometryGrid), intent(in) :: self
    type(coordList), intent(inout) :: coords
    real(defReal), intent(inout)   :: maxDist
    integer(shortInt), intent(out) :: event
    type(distCache), intent(inout) :: cache
    character(100), parameter :: Here = 'move_withCache (geometryGrid_class.f90)'

    ! Unnecessary
    call fatalError(Here, 'Should not be called')

  end subroutine move_withCache

  !!
  !! Move a particle in the top (global) level in the geometry
  !!
  !! See geometry_inter for details
  !!
  !! Uses explicit BC
  !!
  subroutine moveGlobal(self, coords, maxDist, event)
    class(geometryGrid), intent(in) :: self
    type(coordList), intent(inout) :: coords
    real(defReal), intent(inout)   :: maxDist
    integer(shortInt), intent(out) :: event
    real(defReal)                  :: dist
    real(defReal), dimension(3)    :: r, u, r_bar, geomSize

    ! Calculate distance to next cell crossing
    r = coords % lvl(1) % r
    u = coords % lvl(1) % dir
    geomSize = self % geomBounds(4:6) - self % geomBounds(1:3)
    r_bar = -r + self % corner + (HALF + sign(HALF, u)) * geomSize
    dist = minval(r_bar / u)

    if (maxDist < dist) then ! Moves within geometry
      call coords % moveGlobal(maxDist)
      event = COLL_EV
      maxDist = maxDist ! Left for explicitness and compiler

      ! Place coords back into geometry
      call self % placeCoord(coords)

    else ! Hit geometry bounds, increased by NUDGE to avoid numerical issues
      call coords % moveGlobal(dist + NUDGE)
      event = BOUNDARY_EV
      maxDist = dist + NUDGE

      ! Apply boundary conditions 
      call self % explicitBC(coords)

    end if


  end subroutine moveGlobal

  !!
  !! Move a particle in the top level without stopping
  !!
  !! See geometry_inter for details
  !!
  !! Uses co-ordinate transform boundary XSs
  !!
  subroutine teleport(self, coords, dist)
    class(geometryGrid), intent(in) :: self
    type(coordList), intent(inout) :: coords
    real(defReal), intent(in)      :: dist

    ! Move the coords above the geometry
    call coords % moveGlobal(dist)

    ! Place coordinates back into geometry
    call self % placeCoord(coords)

    ! If point is outside apply boundary transformations
    if (coords % matIdx == OUTSIDE_MAT) then
      call self % explicitBC(coords)
    end if

  end subroutine teleport

  !!
  !! CSG equivalent normally in surface class, since we are not using surfaces perform boundary
  !! transformations here instead
  !!
  !! Note indexing difference between bounds and boundary:
  !!   bounds   = [xmin, ymin, zmin, xmax, ...]
  !!   boundary = [x-,   x+,   y-,   y+, z-,z+]
  !!
  subroutine explicitBC(self, coords)
    class(geometryGrid), intent(in) :: self
    class(coordList), intent(inout) :: coords
    integer(shortInt)               :: i
    real(defReal)                   :: outside, move

    ! Loop through axes
    do i = 1, 3
      move = ZERO

      ! Negative side of bounding box
      outside = (self % geomBounds(i) - coords % lvl(1) % r(i))
      if (outside >= ZERO .and. self % boundary(2*i-1) == 1) move = outside

      ! Positive side of bounding box
      outside = (coords % lvl(1) % r(i) - self % geomBounds(i+3)) 
      if (outside >= ZERO .and. self % boundary(2*i) == 1) move = outside

      ! Move if necessary
      if (move > ZERO) then
        ! Flip direction
        coords % lvl(1) % dir(i) = -coords % lvl(1) % dir(i)
        ! Move back into geometry
        call self % teleport(coords, 2*move/abs(coords % lvl(1) % dir(i)))
      end if

    end do

  end subroutine explicitBC

  !!
  !! Returns the list of active materials used in the geometry
  !!
  !! See geometry_inter for details
  !!
  !! NOTE: This function uses VOID_MAT and UNDEF_MAT from universalVariables
  !!       VOID_MAT will be listed multiple times if it occurs in multiple locations
  !!
  function activeMats(self) result(matList)
    class(geometryGrid), intent(in)              :: self
    integer(shortInt), dimension(:), allocatable :: matList

    !TODO: For some reason this gives a warning after compiling, can't figure out why?
    matList = reshape(self % mats,[size(self % mats)])

  end function activeMats

  !!
  !! Return position bounds of a material matIdx
  !!
  !! Result:
  !!   matBounds -> [xmin,ymin,zmin,xmax,ymax,zmax]
  !!
  function matBounds(self, matIdx)
    class(geometryGrid), intent(in) :: self
    integer(shortInt), intent(in)   :: matIdx
    real(defReal), dimension(6)     :: matBounds
    integer(shortInt), dimension(3) :: ijk
    integer(shortInt)               :: i, localID, temp, base
    character(nameLen)              :: matName
    character(100), parameter       :: Here = 'matBounds (geometryGrid_class.f90)'

    ! Convert matIdx to positional ID using name (will be different if void regions exist)
    matName = mm_matName(matIdx)
    read (matName(4:), '(I10)') localID

    ! Convert positional ID to ijk position (same as in latUniverse_class)
    temp = localID - 1
    base = temp / self % latSizeN(1)
    ijk(1) = temp - self % latSizeN(1) * base + 1
    temp = base
    base = temp / self % latSizeN(2)
    ijk(2) = temp - self % latSizeN(2) * base + 1
    ijk(3) = base + 1

    ! Confirm that position is correct
    if (self % mats(ijk(1),ijk(2),ijk(3)) /= matIdx) then
      call fatalError(Here, 'Obtained matIdx different to requested matIdx')
    end if

    ! Set bounds of material
    do i=1, 3
      matBounds(i)   = (ijk(i)-1) * self % latPitch(i) + self % geomBounds(i) + SURF_TOL
      matBounds(i+3) = ijk(i)     * self % latPitch(i) + self % geomBounds(i) - SURF_TOL
    end do

  end function matBounds
    

!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! CSG procedures used for initialisation
!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Descend down the geometry structure untill material is reached
  !!
  !! Requires strting level to be specified.
  !! It is private procedure common to all movment types in geometry.
  !!
  !! Args:
  !!   coords [inout] -> CoordList of a particle. Assume thet coords are already valid for all
  !!     levels above and including start
  !!   start [in] -> Starting level for meterial search
  !!
  !! Errors:
  !!   fatalError if material cell is not found untill maximum nesting is reached
  !!
  subroutine csg_diveToMat(self, coords, start)
    class(geometryGrid), intent(in) :: self
    type(coordList), intent(inout) :: coords
    integer(shortInt), intent(in)  :: start
    integer(shortInt)              :: rootID, localID, fill, id, i
    class(universe), pointer       :: uni
    real(defReal), dimension(3)    :: offset
    character(100), parameter :: Here = 'diveToMat (geometryGrid_class.f90)'

    do i = start, HARDCODED_MAX_NEST
      ! Find cell fill
      rootId = coords % lvl(i) % uniRootID
      localID = coords % lvl(i) % localID
      call self % geom % graph % getFill(fill, id, rootID, localID)

      if (fill >= 0) then ! Found material cell
        coords % matIdx   = fill
        coords % uniqueID = id
        return

      else ! Universe fill descend a level
        if (i == HARDCODED_MAX_NEST) exit ! If there is nested universe at the lowest level

        fill = abs(fill)

        ! Get current universe
        uni => self % geom % unis % getPtr_fast(coords % lvl(i) % uniIdx)

        ! Get cell offset
        offset = uni % cellOffset(coords % lvl(i))

        ! Get nested universe
        uni => self % geom % unis % getPtr_fast(fill)

        ! Enter nested univers
        call coords % addLevel()
        call uni % enter(coords % lvl(i+1), coords % lvl(i) % r - offset, coords % lvl(i) % dir)
        coords % lvl(i+1) % uniRootID = id ! Must be after enter where coord has intent out

      end if
    end do

    call fatalError(Here, 'Failed to find material cell. Should not happen after &
                          &geometry checks during build...')

  end subroutine csg_diveToMat


  !!
  !! Place coordinate list into geometry
  !!
  !! See geometry_inter for details
  !!
  subroutine csg_placeCoord(self, coords)
    class(geometryGrid), intent(in) :: self
    type(coordList), intent(inout) :: coords
    class(universe), pointer       :: uni
    real(defReal), dimension(3)    :: r, dir
    character(100), parameter :: Here = 'placeCoord (geometryGrid_class.f90)'

    ! Check that coordList is initialised
    if (coords % nesting < 1) then
      call fatalError(Here, 'CoordList is not initialised. Nesting is: '//&
                             numToChar(coords % nesting))
    end if

    ! Place coordinates above geometry (in case they were placed)
    call coords % takeAboveGeom()

    ! Enter root universe
    r = coords % lvl(1) % r
    dir = coords % lvl(1) % dir

    uni => self % geom % unis % getPtr_fast(self % geom % rootIdx)

    call uni % enter(coords % lvl(1), r, dir)

    coords % lvl(1) % uniRootID = 1

    ! Dive to material
    call self % csg_diveToMat(coords, 1)

  end subroutine csg_placeCoord


  !!
  !! Find material and unique cell at a given location
  !!
  !! See geometry_inter for details
  !!
  subroutine csg_whatIsAt(self, matIdx, uniqueID, r, u)
    class(geometryGrid), intent(in)                   :: self
    integer(shortInt), intent(out)                    :: matIdx
    integer(shortInt), intent(out)                    :: uniqueID
    real(defReal), dimension(3), intent(in)           :: r
    real(defReal), dimension(3), optional, intent(in) :: u
    type(coordList)                                   :: coords
    real(defReal), dimension(3)                       :: u_l

    ! Select direction
    if (present(u)) then
      u_l = u
    else
      u_l = [ONE, ZERO, ZERO]
    end if

    ! Initialise coordinates
    call coords % init(r, u_l)

    ! Place coordinates
    call self % csg_placeCoord(coords)

    ! Return material & uniqueID
    matIdx   = coords % matIdx
    uniqueID = coords % uniqueID

  end subroutine csg_whatIsAt


end module geometryGrid_class
