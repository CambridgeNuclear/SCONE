module cellUniverse_class

  use numPrecision
  use genericProcedures, only : fatalError, targetNotFound
  use vector_class,      only : vector
  use coord_class,       only : coord
  use surface_inter,     only : surfaceShelf
  use cell_class,        only : cell, cellSHelf
  use universe_inter,    only : universe

  implicit none
  private

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
  subroutine init(self, offset, cellIDs, cShelf)
    class(cellUniverse), intent(inout)        :: self
    real(defReal),dimension(3), intent(in)    :: offset
    integer(shortInt),dimension(:),intent(in) :: cellIDs
    type(cellShelf), intent(in)               :: cShelf
    character(100),parameter :: Here = 'init (cellUniverse_class.f90)'

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

    ! For now there is nothing fancy just try to locate next cell within universe
    call self % findCell(coords, cShelf, sShelf)

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
