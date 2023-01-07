module geomGrid_class

  use numPrecision
  use universalVariables
  use genericProcedures,  only : fatalError
  use dictionary_class,   only : dictionary
  use charMap_class,      only : charMap
  use coord_class,        only : coordList
  use geometry_inter,     only : geometry, distCache
  use csg_class,          only : csg
  use universe_inter,     only : universe

  type, public, extends(geometry) :: geomGrid
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
    procedure :: moveGlobal
    procedure :: teleport
    procedure :: activeMats

    procedure, private :: closestDist

  end type geomGrid

contains

  !!
  !! Initialise geometry
  !!
  !! See geometry_inter for details
  !!
  subroutine init(self, dict, mats, silent)
    class(geomGrid), intent(inout)         :: self
    class(dictionary), intent(in)          :: dict
    type(charMap), intent(in)              :: mats
    logical(defBool), optional, intent(in) :: silent
    logical(defBool)                       :: loud

    ! Build the representation
    call self % geom % init(dict, mats, silent)

  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(geomGrid), intent(inout) :: self

    call self % geom % kill()

  end subroutine kill

  !!
  !! Simply returns distance to closest surface in maxDist variable
  !!
  !! Unlike in geometryStd, here we do not move particle at all
  !!
  subroutine move_noCache(self, coords, maxDist, event)
    class(geomGrid), intent(in)    :: self
    type(coordList), intent(inout) :: coords
    real(defReal), intent(inout)   :: maxDist
    integer(shortInt), intent(out) :: event
    integer(shortInt)              :: surfIdx, lvl
    character(100), parameter :: Here = 'move_noCache (geomGrid_class.f90)'

    call self % closestDist(maxDist, surfIdx, lvl, coords)

  end subroutine move_noCache

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
    class(geomGrid), intent(in)    :: self
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
  !! Unused superclass procedures
  !!

  subroutine placeCoord(self, coords)
    class(geomGrid), intent(in)    :: self
    type(coordList), intent(inout) :: coords
    character(100), parameter :: Here = 'placeCoord (geomGrid_class.f90)'

    call fatalError(Here, "Should not be called")

  end subroutine placeCoord

  subroutine whatIsAt(self, matIdx, uniqueID, r, u)
    class(geomGrid), intent(in)                       :: self
    integer(shortInt), intent(out)                    :: matIdx
    integer(shortInt), intent(out)                    :: uniqueID
    real(defReal), dimension(3), intent(in)           :: r
    real(defReal), dimension(3), optional, intent(in) :: u
    character(100), parameter :: Here = 'whatIsAt (geomGrid_class.f90)'

    call fatalError(Here, "Should not be called")

  end subroutine whatIsAt

  function bounds(self)
    class(geomGrid), intent(in) :: self
    real(defReal), dimension(6) :: bounds
    character(100), parameter :: Here = 'bounds (geomGrid_class.f90)'

    call fatalError(Here, "Should not be called")

  end function bounds

  subroutine move_withCache(self, coords, maxDist, event, cache)
    class(geomGrid), intent(in)    :: self
    type(coordList), intent(inout) :: coords
    real(defReal), intent(inout)   :: maxDist
    integer(shortInt), intent(out) :: event
    type(distCache), intent(inout) :: cache
    character(100), parameter :: Here = 'move_withCache (geomGrid_class.f90)'

    call fatalError(Here, "Should not be called")

  end subroutine move_withCache

  subroutine moveGlobal(self, coords, maxDist, event)
    class(geomGrid), intent(in)    :: self
    type(coordList), intent(inout) :: coords
    real(defReal), intent(inout)   :: maxDist
    integer(shortInt), intent(out) :: event
    character(100), parameter :: Here = 'moveGlobal (geomGrid_class.f90)'

    call fatalError(Here, "Should not be called")

  end subroutine moveGlobal

  subroutine teleport(self, coords, dist)
    class(geomGrid), intent(in)    :: self
    type(coordList), intent(inout) :: coords
    real(defReal), intent(in)      :: dist
    character(100), parameter :: Here = 'teleport (geomGrid_class.f90)'

    call fatalError(Here, "Should not be called")

  end subroutine teleport

  function activeMats(self) result(matList)
    class(geomGrid), intent(in)                  :: self
    integer(shortInt), dimension(:), allocatable :: matList
    character(100), parameter :: Here = 'activeMats (geomGrid_class.f90)'

    call fatalError(Here, "Should not be called")

  end function activeMats

end module geomGrid_class
