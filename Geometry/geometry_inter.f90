module geometry_inter

  use numPrecision
  use universalVariables, only : X_AXIS, Y_AXIS, Z_AXIS, HARDCODED_MAX_NEST
  use genericProcedures,  only : fatalError
  use dictionary_class,   only : dictionary
  use charMap_class,      only : charMap
  use coord_class,        only : coordList

  implicit none
  private

  !!
  !! Distance Cache
  !!
  !! Args:
  !!   lvl  -> Current valid level (e.g. 3 means that levels 1-3 have valid distance)
  !!   dist -> Array of distances
  !!   surf -> Array of coresponding surface mementos for crossing
  !!
  type, public :: distCache
    integer(shortInt) :: lvl = 0
    real(defReal), dimension(HARDCODED_MAX_NEST)     :: dist = ZERO
    integer(shortInt), dimension(HARDCODED_MAX_NEST) :: surf = 0
  end type distCache

  !!
  !! Abstract interface for all geometry implementation
  !!
  type, public, abstract :: geometry
  contains
    ! Generic procedures
    generic :: move  => move_noCache, move_withCache, moveRay_withCache, moveRay_noCache

    ! Deferred procedures
    procedure(init), deferred              :: init
    procedure(kill), deferred              :: kill
    procedure(placeCoord), deferred        :: placeCoord
    procedure(whatIsAt), deferred          :: whatIsAt
    procedure(bounds), deferred            :: bounds
    procedure(move_noCache), deferred      :: move_noCache
    procedure(move_withCache), deferred    :: move_withCache
    procedure(moveRay_noCache), deferred   :: moveRay_noCache
    procedure(moveRay_withCache), deferred :: moveRay_withCache
    procedure(moveGlobal), deferred        :: moveGlobal
    procedure(teleport), deferred          :: teleport
    procedure(activeMats), deferred        :: activeMats
    procedure(numberOfCells), deferred     :: numberOfCells

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
    !!   If geometry is infinite in a given axis direction * then *_min = *_max = ZERO
    !!
    function bounds(self)
      import :: geometry, defReal
      class(geometry), intent(in) :: self
      real(defReal), dimension(6) :: bounds
    end function bounds

    !!
    !! Given coordinates placed in the geometry move point through the geometry
    !!
    !! Move by up to maxDist stopping at domain boundary or untill matIdx or uniqueID changes.
    !! When particle hits boundary, boundary conditions are applied before returning.
    !!
    !! Following events can be returned:
    !!   COLL_EV      -> Particle moved by entire maxDist. Collision happens
    !!   BOUNDARY_EV  -> Particle hit domain boundary
    !!   CROSS_EV     -> Particle crossed to a region with different matIdx or uniqueID
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
    subroutine move_noCache(self, coords, maxDist, event)
      import :: geometry, coordList, defReal, shortInt
      class(geometry), intent(in)    :: self
      type(coordList), intent(inout) :: coords
      real(defReal), intent(inout)   :: maxDist
      integer(shortInt), intent(out) :: event
    end subroutine move_noCache

    !!
    !! Given coordinates placed in the geometry move point through the geometry
    !!
    !! Move by up to maxDist stopping at domain boundary or until matIdx or uniqueID changes.
    !! When particle hits boundary, boundary conditions are applied before returning.
    !!
    !! Use distance cache to avoid needless recalculation of the next crossing at
    !! higher levels.
    !!
    !! Following events can be returned:
    !!   COLL_EV      -> Particle moved by entire maxDist. Collision happens
    !!   BOUNDARY_EV  -> Particle hit domain boundary
    !!   CROSS_EV     -> Particle crossed to a region with different matIdx or uniqueID
    !!   LOST_EV      -> Something gone wrong in tracking and particle is lost
    !!
    !! Args:
    !!   coords [inout]  -> Coordinate list of the particle to be moved through the geometry
    !!   maxDict [inout] -> Maximum distance to move the position. If movement is stopped
    !!     prematurely (e.g. hitting boundary), maxDist is set to the distance the particle has
    !!     moved by.
    !!   event [out] -> Event flag that specifies what finished the movement.
    !!
    !! Errors:
    !!   If coords is not placed in the geometry behaviour is unspecified
    !!   If maxDist < 0.0 behaviour is unspecified
    !!
    subroutine move_withCache(self, coords, maxDist, event, cache)
      import :: geometry, coordList, defReal, shortInt, distCache
      class(geometry), intent(in)    :: self
      type(coordList), intent(inout) :: coords
      real(defReal), intent(inout)   :: maxDist
      integer(shortInt), intent(out) :: event
      type(distCache), intent(inout) :: cache
    end subroutine move_withCache

    !!
    !! Move, but ensuring that vacuum boundaries are treated as reflective and communicating
    !! a vacuum strike back.
    !! This is implemented for handling Random Ray/MoC problems where rays are not terminated
    !! as they strike a vacuum boundary, but are instead reflected.
    !!
    !! Identical in interface to move, with the exception that a flag is included
    !! for identifying a vacuum boundary hit:
    !!
    !! hitVacuum [out] -> false if a vacuum was not hit, true if it was.
    !!
    subroutine moveRay_noCache(self, coords, maxDist, event, hitVacuum)
      import :: geometry, coordList, defReal, shortInt, defBool
      class(geometry), intent(in)    :: self
      type(coordList), intent(inout) :: coords
      real(defReal), intent(inout)   :: maxDist
      integer(shortInt), intent(out) :: event
      logical(defBool), intent(out)  :: hitVacuum
    end subroutine moveRay_noCache

    !!
    !! Move, but ensuring that vacuum boundaries are treated as reflective and communicating
    !! a vacuum strike back.
    !! This is implemented for handling Random Ray/MoC problems where rays are not terminated
    !! as they strike a vacuum boundary, but are instead reflected.
    !!
    !! Identical in interface to move_withCache, with the exception that a flag is included
    !! for identifying a vacuum boundary hit:
    !!
    !! hitVacuum [out] -> false if a vacuum was not hit, true if it was.
    !!
    subroutine moveRay_withCache(self, coords, maxDist, event, cache, hitVacuum)
      import :: geometry, coordList, defReal, shortInt, distCache, defBool
      class(geometry), intent(in)    :: self
      type(coordList), intent(inout) :: coords
      real(defReal), intent(inout)   :: maxDist
      integer(shortInt), intent(out) :: event
      type(distCache), intent(inout) :: cache
      logical(defBool), intent(out)  :: hitVacuum
    end subroutine moveRay_withCache

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
    !!   coords [inout] -> Initialised (but not necessarily placed) coordList for a particle to be
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
    !!   coords [inout] -> Initialised (but not necessarily placed) coordList for a particle to be
    !!   moved. Will become placed on exit.
    !!   dist [in] -> Distance by which to move the particle
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

    !!
    !! Returns the list of active materials used in the geometry
    !!
    !! Args:
    !!   None
    !!
    !! Result:
    !!   Integer list with the IDs of the active materials. Void is not considered
    !!   an active materials, even if it can be present in the geometry.
    !!
    function activeMats(self) result(matList)
      import :: geometry, shortInt
      class(geometry), intent(in)                  :: self
      integer(shortInt), dimension(:), allocatable :: matList
    end function activeMats

    !!
    !! Returns the number of unique cells present in the geometry
    !!
    !! Args:
    !!   None
    !!
    !! Result:
    !!   Integer of the number of unique (material containing) cells in the geometry.
    !!
    function numberOfCells(self) result(n)
      import :: geometry, shortInt
      class(geometry), intent(in) :: self
      integer(shortInt)           :: n
    end function numberOfCells

  end interface

contains

end module geometry_inter
