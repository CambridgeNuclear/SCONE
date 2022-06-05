module geometryStd_class

  use numPrecision
  use universalVariables
  use genericProcedures,  only : fatalError, numToChar
  use coord_class,        only : coordList, coord
  use dictionary_class,   only : dictionary
  use charMap_class,      only : charMap
  use geometry_inter,     only : geometry, distCache
  use csg_class,          only : csg
  use universe_inter,     only : universe
  use surface_inter,      only : surface

  ! Nuclear Data
  use materialMenu_mod,   only : nMat


  implicit none
  private

  !!
  !! Public Pointer Cast
  !!
  public :: geometryStd_CptrCast

  !!
  !! Standard Geometry Model
  !!
  !! Typical geometry of a MC Neutron Transport code composed of multiple nested
  !! universes.
  !!
  !! Boundary conditions in diffrent movement models are handeled:
  !!   move          -> explicitBC
  !!   moveGlobal    -> explicitBC
  !!   moveRay       -> explicitBC with vacuum handled as reflective
  !!   moveRayGlobal -> explicitBC with vacuum handled as reflective
  !!   teleport      -> Co-ordinate transfrom
  !!
  !! Sample Dictionary Input:
  !!   geometry {
  !!     type geometryStd;
  !!     <csg_class definition>
  !!    }
  !!
  !! Public Members:
  !!   geom -> Representation of geometry by csg_class. Contains all surfaces, cells and universe
  !!     as well as geometry graph and info about root uni and boundary surface.
  !!
  !! Interface:
  !!   Geometry Interface
  !!
  type, public, extends(geometry) :: geometryStd
    type(csg) :: geom

  contains
    ! Superclass procedures
    procedure :: init
    procedure :: kill
    procedure :: placeCoord
    procedure :: whatIsAt
    procedure :: bounds
    procedure :: move_noCache
    procedure :: move_withCache
    procedure :: moveRay_noCache
    procedure :: moveRay_withCache
    procedure :: moveGlobal
    procedure :: teleport
    procedure :: activeMats
    procedure :: numberOfCells

    ! Private procedures
    procedure, private :: diveToMat
    procedure, private :: closestDist
    procedure, private :: closestDist_cache
  end type geometryStd

contains

  !!
  !! Initialise geometry
  !!
  !! See geometry_inter for details
  !!
  subroutine init(self, dict, mats, silent)
    class(geometryStd), intent(inout)      :: self
    class(dictionary), intent(in)          :: dict
    type(charMap), intent(in)              :: mats
    logical(defBool), optional, intent(in) :: silent

    ! Build the representation
    call self % geom % init(dict, mats, silent)

  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(geometryStd), intent(inout) :: self

    call self % geom % kill()

  end subroutine kill

  !!
  !! Place coordinate list into geometry
  !!
  !! See geometry_inter for details
  !!
  subroutine placeCoord(self, coords)
    class(geometryStd), intent(in) :: self
    type(coordList), intent(inout) :: coords
    class(universe), pointer       :: uni
    real(defReal), dimension(3)    :: r, dir
    character(100), parameter :: Here = 'placeCoord (geometryStd_class.f90)'

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
    call self % diveToMat(coords, 1)

  end subroutine placeCoord

  !!
  !! Find material and unique cell at a given location
  !!
  !! See geometry_inter for details
  !!
  subroutine whatIsAt(self, matIdx, uniqueID, r, u)
    class(geometryStd), intent(in)                    :: self
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
    call self % placeCoord(coords)

    ! Return material & uniqueID
    matIdx   = coords % matIdx
    uniqueID = coords % uniqueID

  end subroutine whatIsAt

  !!
  !! Return Axis Aligned Bounding Box encompassing the geometry
  !!
  !! See geometry_inter for details
  !!
  function bounds(self)
    class(geometryStd), intent(in) :: self
    real(defReal), dimension(6) :: bounds
    class(surface), pointer     :: surf
    integer(shortInt)           :: i

    ! Get boundary surface
    surf => self % geom % surfs % getPtr(self % geom % borderIdx)
    bounds = surf % boundingBox()

    ! Change the infinate dimension to 0.0
    do i = 1, 3
      if(bounds(i) <= -INF .and. bounds(i+3) >= INF) then
        bounds(i) = ZERO
        bounds(i+3) = ZERO
      end if
    end do

  end function bounds

  !!
  !! Given coordinates placed in the geometry move point through the geometry
  !!
  !! See geometry_inter for details
  !!
  !! Uses explicit BC
  !!
  subroutine move_noCache(self, coords, maxDist, event)
    class(geometryStd), intent(in) :: self
    type(coordList), intent(inout) :: coords
    real(defReal), intent(inout)   :: maxDist
    integer(shortInt), intent(out) :: event
    integer(shortInt)              :: surfIdx, level
    real(defReal)                  :: dist
    class(surface), pointer        :: surf
    class(universe), pointer       :: uni
    character(100), parameter :: Here = 'move (geometryStd_class.f90)'

    if (.not.coords % isPlaced()) then
      call fatalError(Here, 'Coordinate list is not placed in the geometry')
    end if

    ! Find distance to the next surface
    call self % closestDist(dist, surfIdx, level, coords)

    if (maxDist < dist) then ! Moves within cell
      call coords % moveLocal(maxDist, coords % nesting)
      event = COLL_EV
      maxDist = maxDist ! Left for explicitness. Compiler will not stand it anyway

    else if (surfIdx == self % geom % borderIdx .and. level == 1) then ! Hits domain boundary
      ! Move global to the boundary
      call coords % moveGlobal(dist)
      event = BOUNDARY_EV
      maxDist = dist

      ! Get boundary surface and apply BCs
      surf => self % geom % surfs % getPtr(self % geom % borderIdx)
      call surf % explicitBC(coords % lvl(1) % r, coords % lvl(1) % dir)

      ! Place back in geometry
      call self % placeCoord(coords)

    else ! Crosses to diffrent local cell
      ! Move to boundary at hit level
      call coords % moveLocal(dist, level)
      event = CROSS_EV
      maxDist = dist

      ! Get universe and cross to the next cell
      uni => self % geom % unis % getPtr_fast(coords % lvl(level) % uniIdx)
      call uni % cross(coords % lvl(level), surfIdx)

      ! Get material
      call self % diveToMat(coords, level)

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
    class(geometryStd), intent(in) :: self
    type(coordList), intent(inout) :: coords
    real(defReal), intent(inout)   :: maxDist
    integer(shortInt), intent(out) :: event
    type(distCache), intent(inout) :: cache
    integer(shortInt)              :: surfIdx, level
    real(defReal)                  :: dist
    class(surface), pointer        :: surf
    class(universe), pointer       :: uni
    character(100), parameter :: Here = 'move_withCache (geometryStd_class.f90)'

    if (.not.coords % isPlaced()) then
      print *, coords % lvl(1) % r
      print *, coords % lvl(1) % dir
      call fatalError(Here, 'Coordinate list is not placed in the geometry')
    end if

    ! Find distance to the next surface
    call self % closestDist_cache(dist, surfIdx, level, coords, cache)

    if (maxDist < dist) then ! Moves within cell
      call coords % moveLocal(maxDist, coords % nesting)
      event = COLL_EV
      maxDist = maxDist ! Left for explicitness. Compiler will not stand it anyway
      cache % lvl = 0

    else if (surfIdx == self % geom % borderIdx .and. level == 1) then ! Hits domain boundary
      ! Move global to the boundary
      call coords % moveGlobal(dist)
      event = BOUNDARY_EV
      maxDist = dist
      cache % lvl = 0

      ! Get boundary surface and apply BCs
      surf => self % geom % surfs % getPtr(self % geom % borderIdx)
      call surf % explicitBC(coords % lvl(1) % r, coords % lvl(1) % dir)

      ! Place back in geometry
      call self % placeCoord(coords)

    else ! Crosses to diffrent local cell
      ! Move to boundary at hit level
      call coords % moveLocal(dist, level)
      event = CROSS_EV
      maxDist = dist
      cache % dist(1:level-1) = cache % dist(1:level-1) - dist
      cache % lvl = level - 1

      ! Get universe and cross to the next cell
      uni => self % geom % unis % getPtr_fast(coords % lvl(level) % uniIdx)
      call uni % cross(coords % lvl(level), surfIdx)

      ! Get material
      call self % diveToMat(coords, level)

    end if

  end subroutine move_withCache
  
  !!
  !! Given coordinates placed in the geometry move point through the geometry
  !!
  !! See geometry_inter for details
  !!
  !! Uses explicit BC while treating vacuum boundaries as reflective
  !!
  subroutine moveRay_noCache(self, coords, maxDist, event, hitVacuum)
    class(geometryStd), intent(in) :: self
    type(coordList), intent(inout) :: coords
    real(defReal), intent(inout)   :: maxDist
    integer(shortInt), intent(out) :: event
    logical(defBool), intent(out)  :: hitVacuum
    integer(shortInt)              :: surfIdx, level
    real(defReal)                  :: dist
    class(surface), pointer        :: surf
    class(universe), pointer       :: uni
    character(100), parameter :: Here = 'moveRay_noCache (geometryStd_class.f90)'

    hitVacuum = .FALSE.

    if (.not.coords % isPlaced()) then
      print *, coords % lvl(1) % r
      print *, coords % lvl(1) % dir
      call fatalError(Here, 'Coordinate list is not placed in the geometry')
    end if

    ! Find distance to the next surface
    call self % closestDist(dist, surfIdx, level, coords)

    if (maxDist < dist) then ! Moves within cell
      call coords % moveLocal(maxDist, coords % nesting)
      event = COLL_EV
      maxDist = maxDist ! Left for explicitness. Compiler will not stand it anyway

    else if (surfIdx == self % geom % borderIdx .and. level == 1) then ! Hits domain boundary
      ! Move global to the boundary
      call coords % moveGlobal(dist)
      event = BOUNDARY_EV
      maxDist = dist

      ! Get boundary surface and apply BCs
      surf => self % geom % surfs % getPtr(self % geom % borderIdx)
      call surf % explicitRayBC(coords % lvl(1) % r, coords % lvl(1) % dir, hitVacuum)

      ! Place back in geometry
      call self % placeCoord(coords)

    else ! Crosses to diffrent local cell
      ! Move to boundary at hit level
      call coords % moveLocal(dist, level)
      event = CROSS_EV
      maxDist = dist

      ! Get universe and cross to the next cell
      uni => self % geom % unis % getPtr_fast(coords % lvl(level) % uniIdx)
      call uni % cross(coords % lvl(level), surfIdx)

      ! Get material
      call self % diveToMat(coords, level)

    end if

  end subroutine moveRay_noCache


  !!
  !! Given coordinates placed in the geometry move point through the geometry
  !!
  !! See geometry_inter for details
  !!
  !! Uses explicit BC while treating vacuum boundaries as reflective
  !!
  subroutine moveRay_withCache(self, coords, maxDist, event, cache, hitVacuum)
    class(geometryStd), intent(in) :: self
    type(coordList), intent(inout) :: coords
    real(defReal), intent(inout)   :: maxDist
    integer(shortInt), intent(out) :: event
    type(distCache), intent(inout) :: cache
    logical(defBool), intent(out)  :: hitVacuum
    integer(shortInt)              :: surfIdx, level
    real(defReal)                  :: dist
    class(surface), pointer        :: surf
    class(universe), pointer       :: uni
    character(100), parameter :: Here = 'moveRay_withCache (geometryStd_class.f90)'

    hitVacuum = .FALSE.

    if (.not.coords % isPlaced()) then
      print *, coords % lvl(1) % r
      print *, coords % lvl(1) % dir
      call fatalError(Here, 'Coordinate list is not placed in the geometry')
    end if

    ! Find distance to the next surface
    call self % closestDist_cache(dist, surfIdx, level, coords, cache)

    if (maxDist < dist) then ! Moves within cell
      call coords % moveLocal(maxDist, coords % nesting)
      event = COLL_EV
      maxDist = maxDist ! Left for explicitness. Compiler will not stand it anyway
      cache % lvl = 0

    else if (surfIdx == self % geom % borderIdx .and. level == 1) then ! Hits domain boundary
      ! Move global to the boundary
      call coords % moveGlobal(dist)
      event = BOUNDARY_EV
      maxDist = dist
      cache % lvl = 0

      ! Get boundary surface and apply BCs
      surf => self % geom % surfs % getPtr(self % geom % borderIdx)
      call surf % explicitRayBC(coords % lvl(1) % r, coords % lvl(1) % dir, hitVacuum)

      ! Place back in geometry
      call self % placeCoord(coords)

    else ! Crosses to diffrent local cell
      ! Move to boundary at hit level
      call coords % moveLocal(dist, level)
      event = CROSS_EV
      maxDist = dist
      cache % dist(1:level-1) = cache % dist(1:level-1) - dist
      cache % lvl = level - 1

      ! Get universe and cross to the next cell
      uni => self % geom % unis % getPtr_fast(coords % lvl(level) % uniIdx)
      call uni % cross(coords % lvl(level), surfIdx)

      ! Get material
      call self % diveToMat(coords, level)

    end if

  end subroutine moveRay_withCache

  !!
  !! Move a particle in the top (global) level in the geometry
  !!
  !! See geometry_inter for details
  !!
  !! Uses explicit BC
  !!
  subroutine moveGlobal(self, coords, maxDist, event)
    class(geometryStd), intent(in) :: self
    type(coordList), intent(inout) :: coords
    real(defReal), intent(inout)   :: maxDist
    integer(shortInt), intent(out) :: event
    class(surface), pointer        :: surf
    real(defReal)                  :: dist

    ! Get boundary surface
    surf => self % geom % surfs % getPtr(self % geom % borderIdx)

    ! Find distance to the boundary
    dist = surf % distance(coords % lvl(1) % r, coords % lvl(1) % dir)

    ! Select collision or boundary hit
    if (maxDist < dist) then ! maxDist is shorter
      ! Move above the geometry
      call coords % moveGlobal(maxDist)
      event = COLL_EV
      maxDist = maxDist

    else
      ! Move to boundary and apply BC
      call coords % moveGlobal(dist)
      event = BOUNDARY_EV
      call surf % explicitBC(coords % lvl(1) % r, coords % lvl(1) % dir)
      maxDist = dist

    end if

    ! Return particle to geometry
    call self % placeCoord(coords)

  end subroutine moveGlobal

  !!
  !! Move a particle in the top level without stopping
  !!
  !! See geometry_inter for details
  !!
  !! Uses co-ordinate transform boundary XSs
  !!
  subroutine teleport(self, coords, dist)
    class(geometryStd), intent(in) :: self
    type(coordList), intent(inout) :: coords
    real(defReal), intent(in)      :: dist
    class(surface), pointer        :: surf

    ! Move the coords above the geometry
    call coords % moveGlobal(dist)

    ! Place coordinates back into geometry
    call self % placeCoord(coords)

    ! If point is outside apply boundary transformations
    if (coords % matIdx == OUTSIDE_MAT) then
      surf => self % geom % surfs % getPtr(self % geom % borderIdx)
      call surf % transformBC(coords % lvl(1) % r, &
                              coords % lvl(1) % dir)

      ! Return particle to geometry
      call self % placeCoord(coords)
    end if

  end subroutine teleport

  !!
  !! Returns the list of active materials used in the geometry
  !!
  !! See geometry_inter for details
  !!
  !! NOTE: This function uses VOID_MAT and UNDEF_MAT from universalVariables
  !!
  function activeMats(self) result(matList)
    class(geometryStd), intent(in)               :: self
    integer(shortInt), dimension(:), allocatable :: matList
    integer(shortInt)                            :: N, lastIdx

    ! Takes the list of materials present in the geometry from geomGraph
    N = size(self % geom % graph % usedMats)
    lastIdx = self % geom % graph % usedMats(N)

    ! Check if the last entry of the list is an actual material or void
    if (lastIdx == VOID_MAT) then
      N = N - 1
      lastIdx = self % geom % graph % usedMats(N)
    end if
    ! Check if the last entry of the list is an undefined material
    if (lastIdx == UNDEF_MAT) then
      matList = self % geom % graph % usedMats(1:N-1)
    else
      matList = self % geom % graph % usedMats(1:N)
    end if

  end function activeMats
  
  !!
  !! Returns the number of unique cells present in the geometry
  !!
  !! See geometry_inter for details
  !!
  function numberOfCells(self) result(n)
    class(geometryStd), intent(in) :: self
    integer(shortInt)              :: n

    ! Takes the number of cells from geomGraph
    n = self % geom % graph % uniqueCells

  end function numberOfCells

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
  subroutine diveToMat(self, coords, start)
    class(geometryStd), intent(in) :: self
    type(coordList), intent(inout) :: coords
    integer(shortInt), intent(in)  :: start
    integer(shortInt)              :: rootID, localID, fill, id, i
    class(universe), pointer       :: uni
    real(defReal), dimension(3)    :: offset
    character(100), parameter :: Here = 'diveToMat (geometryStd_class.f90)'

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

  end subroutine diveToMat

  !!
  !! Return distance to the closest surface
  !!
  !! Searches through all geometry levels. In addition to distance return level
  !! and surfIdx for crossing surface
  !!
  !! Args:
  !!   dist [out]    -> Value of closest distance
  !!   surfIdx [out] -> Surface index for the crossing returned from the universe
  !!   lvl     [out] -> Level at which crossing is closest
  !!   coords [in]   -> Current coordinates of a particle
  !!
  subroutine closestDist(self, dist, surfIdx, lvl, coords)
    class(geometryStd), intent(in) :: self
    real(defReal), intent(out)     :: dist
    integer(shortInt), intent(out) :: surfIdx
    integer(shortInt), intent(out) :: lvl
    type(coordList), intent(in)    :: coords
    integer(shortInt)              :: l, test_idx
    real(defReal)                  :: test_dist
    class(universe), pointer       :: uni

    dist = INF
    surfIdx = 0
    lvl = 0
    do l = 1, coords % nesting
      ! Get universe
      uni => self % geom % unis % getPtr_fast(coords % lvl(l) % uniIdx)

      ! Find distance
      call uni % distance(test_dist, test_idx, coords % lvl(l))

      ! Save distance, surfIdx & level coresponding to shortest distance
      ! Take FP precision into account
      if ((dist - test_dist) >= dist * FP_REL_TOL) then
        dist = test_dist
        surfIdx = test_idx
        lvl = l
      end if

    end do
  end subroutine closestDist

  !!
  !! Return distance to the closest surface
  !!
  !! Searches through all geometry levels. In addition to distance return level
  !! and surfIdx for crossing surface
  !!
  !! Args:
  !!   dist [out]    -> Value of closest distance
  !!   surfIdx [out] -> Surface index for the crossing returned from the universe
  !!   lvl     [out] -> Level at which crossing is closest
  !!   coords [in]   -> Current coordinates of a particle
  !!   cache [inout] -> Distance cache. Use valid distances from cache. Put calculated
  !!     distances on the cache.
  !!
  subroutine closestDist_cache(self, dist, surfIdx, lvl, coords, cache)
    class(geometryStd), intent(in) :: self
    real(defReal), intent(out)     :: dist
    integer(shortInt), intent(out) :: surfIdx
    integer(shortInt), intent(out) :: lvl
    type(coordList), intent(in)    :: coords
    type(distCache), intent(inout) :: cache
    integer(shortInt)              :: l, test_idx
    real(defReal)                  :: test_dist
    class(universe), pointer       :: uni

    dist = INF
    surfIdx = 0
    lvl = 0
    do l = 1, coords % nesting

      ! Update Cache if distance is not valid
      if (cache % lvl < l) then
        ! Get universe
        uni => self % geom % unis % getPtr_fast(coords % lvl(l) % uniIdx)

        ! Find distance
        call uni % distance(cache % dist(l), cache % surf(l), coords % lvl(l))
        cache % lvl = cache % lvl + 1
      end if

      ! Read distance and crossing memento from cache
      test_dist = cache % dist(l)
      test_idx  = cache % surf(l)

      ! Save distance, surfIdx & level coresponding to shortest distance
      ! Take FP precision into account
      if ((dist - test_dist) >= dist * FP_REL_TOL) then
        dist = test_dist
        surfIdx = test_idx
        lvl = l
      end if

    end do
  end subroutine closestDist_cache

  !!
  !! Cast geometry pointer to geometryStd class pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class geometry
  !!
  !! Result:
  !!   Null if source is not of geometryStd class
  !!   Target points to source if source is geometryStd class
  !!
  pure function geometryStd_CptrCast(source) result(ptr)
    class(geometry), pointer, intent(in) :: source
    class(geometryStd), pointer          :: ptr

    select type(source)
      class is(geometryStd)
        ptr => source

      class default
        ptr => null()
    end select

  end function geometryStd_CptrCast

end module geometryStd_class
