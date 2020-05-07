module cellUniverse_class

  use numPrecision
  use universalVariables
  use genericProcedures, only : fatalError, targetNotFound, numToChar
  use vector_class,      only : vector
  use dictionary_class,  only : dictionary
  use coord_class,       only : coord
  use intMap_class,      only : intMap
  use charMap_class,     only : charMap
  use surface_inter,     only : surfaceSlot, surfaceShelf
  use cell_class,        only : cell, cellShelf
  use universe_inter,    only : universe

  implicit none
  private

  !!
  !! Constructor
  !!
  interface cellUniverse
    module procedure cellUniverse_fromDict
  end interface

  !!
  !! Basic implementation of the universe
  !! It is composed only from cells in cellShelf
  !!
  type, public, extends(universe) :: cellUniverse
    private
    integer(shortInt),dimension(:), allocatable :: cellIDXs

  contains
    ! Build procedures
    procedure :: init
    procedure :: cellIdx

    ! Runtime procedures
    procedure :: findCell
    procedure :: distance
    procedure :: cross
    procedure :: cellOffset

  end type cellUniverse

contains

  !!
  !! Initialise cellUniverse
  !!
  subroutine init(self, offset, id, cellIDs, cShelf)
    class(cellUniverse), intent(inout)        :: self
    real(defReal),dimension(3), intent(in)    :: offset
    integer(shortInt), intent(in)             :: id
    integer(shortInt),dimension(:),intent(in) :: cellIDs
    type(cellShelf), intent(in)               :: cShelf
    character(100),parameter :: Here = 'init (cellUniverse_class.f90)'

    ! Load ID
    if( id < 1) call fatalError(Here,'Invalid id: ' // numToChar(id))
    call self % setId(Id)

    ! Load universe offset
    call self% setOffset(offset)

    ! Load cellIdxs
    self % cellIDXs = cShelf % getIdx(cellIDs)

    ! Verify that all cells were present
    if(any(self % cellIDXs == targetNotFound)) then
      call fatalError(Here, 'One of the cells is not defined in cellShelf')
    end if

  end subroutine init

  !!
  !! Returns an initialised instance of cell universe from dict and cellShelf
  !! Return fillVector as well. +ve entries are materials IDXs, -ve are fill universes IDs
  !!
  function cellUniverse_fromDict(fillVector, dict, cShelf, sShelf, cellFillMap, materials) result (new)
    integer(shortInt),dimension(:),allocatable,intent(out) :: fillVector
    class(dictionary), intent(in)                          :: dict
    type(cellShelf), intent(inout)                         :: cShelf
    type(surfaceShelf), intent(inout)                      :: sShelf
    type(intMap), intent(in)                               :: cellFillMap
    type(charMap), intent(in)                              :: materials
    type(cellUniverse)                                     :: new
    real(defReal), dimension(:),allocatable                :: offset
    integer(shortInt)                                      :: id, i
    integer(shortInt),dimension(:), allocatable            :: cellIDs
    character(100), parameter :: Here = ' cellUniverse_fromDict (cellUniverse_class.f90)'

    ! Load all data
    call dict % get(id,'id')
    call dict % get( cellIds,'cells')

    call dict % getOrDefault( offset, 'origin', [ZERO, ZERO, ZERO] )

    ! Verify data
    if( size(offset)  /= 3) call fatalError(Here,'Origin needs to be size 3')
    if( size(cellIds) == 0) call fatalError(Here,'Universe needs to contain some cells')

    ! Initialise universe
    call new % init(offset,id,cellIDs, cShelf)

    ! Allocate fillVector
    fillVector = cellIDs

    ! Translate cell fills to matarial and uniIds
    do i= 1, size(cellIds)
      fillVector(i) = cellFillMap % get(cellIDs(i))

    end do

  end function cellUniverse_fromDict

  !!
  !! Give cellIdx given localID of cell
  !!
  function cellIdx(self, localID)
    class(cellUniverse), intent(in) :: self
    integer(shortInt), intent(in) :: localID
    integer(shortInt)             :: cellIdx

    cellIdx = self % cellIDXs(localID)

  end function cellIdx


  !!
  !! Using the coordinates it finds a localID & cellIDx inside the universe
  !!
  subroutine findCell(self, coords, cShelf, sShelf)
    class(cellUniverse), intent(in)  :: self
    type(coord),intent(inout)        :: coords
    type(cellShelf), intent(in)      :: cShelf
    type(surfaceShelf), intent(in)   :: sShelf
    integer(shortInt)                :: i
    type(vector)                     :: r, u
    character(100), parameter       :: Here = 'findCell (cellUniverse_class.f90)'

    ! Copy position and direction into vectors
    r = coords % r
    u = coords % dir

    ! Loop through all cells
    do i = 1, size(self % cellIDXs)

      ! Associate i-th cell with word cell
      associate ( ith_cell => cShelf % shelf ( self % cellIDXs(i)) )

        ! If point is inside the cell return current cell
        if( ith_cell % isInside(r, u, sShelf) ) then
          coords % cellIdx = self % cellIDXs(i)
          coords % localID = i
          return

        end if
      end associate
    end do

    ! Cell was not found
    call fatalError(Here,'Cell containing point and direction was not found.')

  end subroutine findCell

  !!
  !! Returns distance to the next cell boundary in the universe
  !! Returns surfIdx of the surface beeing X-ed
  !! surfIdx can be < 0 (invalid)
  !! Invalid surfIdx indicate a surface that is private to the universe(not present on sShelf)
  !! Example of a private surface is face of a lattice cell
  !!
  subroutine distance(self, dist, surfIdx, coords ,cShelf, sShelf)
    class(cellUniverse), intent(in)  :: self
    real(defReal), intent(out)       :: dist
    integer(shortInt), intent(out)   :: surfIdx
    type(coord), intent(in)          :: coords
    type(cellShelf), intent(in)      :: cShelf
    type(surfaceShelf),intent(in)    :: sShelf
    integer(shortInt)                :: cellIdx, localID
    type(vector)                     :: r, u
    character(100),parameter         :: Here ='distance (cellUniverse_class.f90)'

    ! Get current cell from the coords
    localID = coords % localID
    cellIdx = coords % cellIdx

    ! Get position and direction
    r = coords % r
    u = coords % dir

    ! Verify that cell Idx and localID agree
    if (cellIdx /= self % cellIDXs(localID) ) then
      call fatalError(Here,'localID and cellIdx on coord clash with universe cellIDXs vector')
    end if

    ! Obtain distance
    call cShelf % shelf(cellIdx) % distance(dist, surfIdx, r, u, sShelf)

  end subroutine distance

  !!
  !! Perform crossing inside the universe from current cell to next cell
  !! Assumes coords are at the surface being crossed (within surface_tol)
  !! IT DOES NOT PERFORM CHECK IF THE ABOVE ASSUMPTION IS VALID! (Be careful - MAK)
  !!
  subroutine cross(self, coords, surfIdx, cShelf, sShelf)
    class(cellUniverse), intent(in) :: self
    type(coord),intent(inout)       :: coords
    integer(shortInt), intent(in)   :: surfIdx
    type(cellShelf), intent(in)     :: cShelf
    type(surfaceShelf), intent(in)  :: sShelf
    real(defReal), dimension(3)     :: r_in

    ! Nudge position slightly forward to avoid dot product with normal vector at the surface
    ! This is done mainly for optimisation. Avoids calculation normal vector.
    r_in = coords % r
    coords % r = r_in + coords % dir * NUDGE

    ! For now there is nothing fancy just try to locate next cell within universe
    call self % findCell(coords, cShelf, sShelf)

    ! Return to original position
    coords % r = r_in

  end subroutine cross

  !!
  !! Return offset for the current cell
  !! This is used when going into nested universe
  !! Total offset of the nested universe is :
  !!  cellOffset + nestedUniverse % offset
  !!
  !! Cell offset needs to be applied before envoing "enter" on the nested universe
  !!
  function cellOffset(self,coords) result(offset)
    class(cellUniverse), intent(in) :: self
    type(coord),intent(inout)       :: coords
    real(defReal),dimension(3)      :: offset

    offset = [0.0, 0.0, 0.0]

  end function cellOffset

end module cellUniverse_class
