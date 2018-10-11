module pinUniverse_class

  use numPrecision
  use universalVariables
  use genericProcedures,   only : fatalError, hasDuplicates, linFind, targetNotFound, swap, numToChar
  use vector_class,        only : vector
  use dictionary_class,    only : dictionary
  use intMap_class,        only : intMap

  use coord_class,         only : coord
  use surface_inter,       only : surface, surfaceSlot, surfaceShelf
  use surfaceFactory_func, only : new_surface
  use cell_class,          only : cell, cellShelf
  use universe_inter,      only : universe

  !*** STAYS HERE ONLY PROVISIONALLY
  use nuclearData_inter, only : nuclearData

  implicit none
  private

  !!
  !! Constructor
  !!
  interface pinUniverse
    module procedure pinUniverse_fromDict
  end interface

  !!
  !! Small private type to group all data related to given annulus
  !!
  type, private :: annulus
    integer(shortInt) :: cellIdx
    integer(shortInt) :: surfIdx
  end type annulus

  !!
  !!
  !!
  type, public,extends(universe) :: pinUniverse
    private
    real(defReal), dimension(:), allocatable :: r
    real(defReal), dimension(:), allocatable :: r2
    type(annulus), dimension(:), allocatable :: data


  contains
    ! Build procedures
    procedure :: init
    procedure :: cellIdx

    ! Runtime procedures
    procedure :: findCell
    procedure :: distance
    procedure :: cross
    procedure :: cellOffset

    ! Private procedures
    procedure, private :: buildAnnuli
    procedure, private :: makeFillVector
  end type pinUniverse

contains

  !!
  !! Initialise pinUniverse
  !! use u<####> syntax for the universe fill
  !! radius = 0.0 specifies outermost annulus !
  !!
  subroutine init(self, fillVector, id, matNames, radii, sShelf, cShelf, materials)
    integer(shortInt),dimension(:),allocatable,intent(out) :: fillVector
    class(pinUniverse), intent(inout)                      :: self
    integer(shortInt), intent(in)                          :: id
    character(nameLen),dimension(:),intent(in)             :: matNames
    real(defReal), dimension(:), intent(in)                :: radii
    type(surfaceShelf), intent(inout)                      :: sShelf
    type(cellShelf),intent(inout)                          :: cShelf
    class(nuclearData), intent(in)                         :: materials
    character(nameLen),dimension(size(radii))              :: names
    integer(shortInt)                                      :: N, i, idx
    character(100), parameter :: Here = 'init (pinUniverse_class.f90)'

    ! Set id
    call self % setId(id)

    ! Check that radii are -ve and that matNames and radii are the same size
    if (any (radii < ZERO)) call fatalError(Here,'-ve radii cannot be provided')
    if (size(radii) /= size(matNames)) call fatalError(Here,'Size of matNames & radii does not match')
    if (hasDuplicates(radii)) call fatalError(Here,'Duplicate values in annuli radii')

    ! Load number of annuli
    N = size(radii)

    ! Allocate storage space and set offset to 0
    call self % setOffset([ZERO, ZERO, ZERO])
    self % r = radii
    allocate(self % data(N))

    ! Copy material names to local copy
    names = matNames

    ! Radii may not be in right order. Sort them together with corresponding material names

    ! Start by processing the outermost element
    idx = linFind(self % r, ZERO)
    if(idx == targetNotFound) call fatalError(Here,'Outermost element with radius 0.0 was not found')

    call swap( self % r(idx),    self % r(N))
    call swap( names(idx), names(N))

    ! Put the rest in the ascending order -> Selection sort for simplicity.
    ! This is not performace critical (not yet at least)
    do i = N-1,1,-1
      idx = maxloc( self % r(1:i), 1 )
      call swap (self % r(idx), self % r(i))
      call swap (names(idx), names(i) )
    end do

    ! Build surfaces and cells
    call self % buildAnnuli(sShelf, cShelf)

    ! Create build vector
    call self % makeFillVector(fillVector, matNames, materials)

    ! Set outermost radius beyond INFINITY
    self % r(size(self % r)) = 1.1 * INFINITY

    ! Store r2
    self % r2 = self % r * self % r

  end subroutine init

  !!
  !! Returns an initialised instance of a pin universe
  !! Returns fillVector as well with +ve entries being materialIDXs and -ve fill universes IDs
  !! Provisionally provide nuclearData *** REPLACE WITH CHAR-INT MAP LATER FOR DECOUPLING
  !!
  function pinUniverse_fromDict(fillVector, dict, cShelf, sShelf, cellFillMap, materials) result (new)
    integer(shortInt),dimension(:),allocatable,intent(out) :: fillVector
    class(dictionary), intent(in)                          :: dict
    type(cellShelf), intent(inout)                         :: cShelf
    type(surfaceShelf), intent(inout)                      :: sShelf
    type(intMap), intent(in)                               :: cellFillMap
    class(nuclearData), intent(in)                         :: materials
    type(pinUniverse)                                      :: new
    integer(shortInt)                                      :: id
    character(nameLen),dimension(:),allocatable            :: keys
    real(defReal),dimension(:),allocatable                 :: radii

    ! Read universe id and data
    call dict % get(id,'id')
    call dict % get(keys,'fills')
    call dict % get(radii,'radii')

    ! Initialise universe
    call new % init(fillVector, id, keys, radii, sShelf, cShelf, materials)

  end function pinUniverse_fromDict

  !!
  !! Return cell index given localID of a cell
  !!
  function cellIdx(self,localId)
    class(pinUniverse), intent(in) :: self
    integer(shortInt), intent(in)  :: localID
    integer(shortInt)              :: cellIdx
    character(100), parameter :: Here = 'cellIdx (pinUniverse_class.f90)'

    if (localID < 1 .or. localID > size(self % r)) then
      call fatalError(Here,'Provided local ID ' // numToChar(localID) // ' is Not present')
    end if

    cellIdx = self % data(localID) % cellIdx

  end function cellIdx

  !!
  !! Using the coordinates it find a localID & cellIDx inside the universe
  !! NOTE: ALONG WITH GENERAL CYLINDER THIS BREAKS ASSAMPTIONS ABOUT SURFACE TOLERANCE
  !!
  subroutine findCell(self, coords, cShelf, sShelf)
    class(pinUniverse), intent(in) :: self
    type(coord), intent(inout)     :: coords
    type(cellShelf), intent(in)    :: cShelf
    type(surfaceSHelf), intent(in) :: sShelf
    real(defReal)                  :: r2
    integer(shortInt)              :: idx, i

    r2 = coords % r(1) * coords % r(1) + coords % r(2) * coords % r(2)

    ! Do a counting search to find aproperiate ring
    ! Avoids branching
    idx = 1
    do i=2,size(self %r)
      if ( r2 > self % r2(i-1)) idx = idx +1

    end do

    coords % localID = idx
    coords % cellIdx = self % data(idx) % cellIdx

  end subroutine findCell

  !!
  !! Returns distance to the next cell boundary in the universe
  !! Returns -1 if crossing is outward, returns -2 if crossing inward
  !!
  subroutine distance(self, dist, surfIdx, coords, cShelf, sShelf)
    class(pinUniverse), intent(in) :: self
    real(defReal), intent(out)     :: dist
    integer(shortInt), intent(out) :: surfIdx
    type(coord), intent(in)        :: coords
    type(cellShelf), intent(in)    :: cShelf
    type(surfaceShelf),intent(in)  :: sShelf
    integer(shortInt)              :: cellIdx, localID
    type(vector)                   :: r, u
    character(100),parameter :: Here = 'distance (pinUniverse_class.f90)'

    ! Get current cell from the coords
    localID = coords % localID
    cellIdx = coords % cellIdx

    ! Get position and direction
    r = coords % r
    u = coords % dir

    ! Verify that cell Idx and localID agree
    if (cellIdx /= self % data(localID) % cellIdx ) then
      call fatalError(Here,'localID and cellIdx on coord clash with universe cellIDXs vector')
    end if

    ! Obtain distance
    call cShelf % shelf(cellIdx) % distance(dist, surfIdx, r, u, sShelf)

    ! Change surfIdx to -1 or -2
    if( surfIdx == self % data(localID) % surfIdx) then
      surfIdx = -1
    else
      surfIdx = -2
    end if

  end subroutine distance

  !!
  !! Perform crossing inside the universe from current cell to next cell
  !! Assumes coords are at the surface being crossed (within surface_tol)
  !! IT DOES NOT PERFORM CHECK IF THE ABOVE ASSUMPTION IS VALID! (Be careful - MAK)
  !!
  subroutine cross(self, coords, surfIdx, cShelf, sShelf)
    class(pinUniverse), intent(in) :: self
    type(coord),intent(inout)      :: coords
    integer(shortInt), intent(in)  :: surfIdx
    type(cellShelf), intent(in)    :: cShelf
    type(surfaceShelf), intent(in) :: sShelf
    character(100),parameter :: Here = 'cross (pinUniverse_class.f90)'

    ! Find next local cell
    select case(surfIdx)
      case(-1)
        coords % localID = coords % localID + 1

      case(-2)
        coords % localID = coords % localID - 1

      case default
        call fatalError(Here,'Unknown surface index')

    end select

    ! Load cellIdx
    coords % cellIdx = self % data( coords % localID) % cellIdx

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
    class(pinUniverse), intent(in) :: self
    type(coord),intent(inout)      :: coords
    real(defReal),dimension(3)     :: offset

    offset = [0.0, 0.0, 0.0]

  end function cellOffset

  !!
  !! Build all surfaces and cells that define pinUniverse
  !!
  subroutine buildAnnuli(self, sShelf, cShelf)
    class(pinUniverse), intent(inout)           :: self
    type(surfaceShelf), intent(inout)           :: sShelf
    type(cellShelf), intent(inout)              :: cShelf
    type(dictionary)                            :: tempDict
    class(surface),allocatable                  :: tempSurf
    type(cell)                                  :: tempCell
    integer(shortInt)                           :: cellId, i
    integer(shortInt),dimension(size(self % r)) :: surfIDs

    ! Build all surfaces
    do i =1,size(self % r)
      ! Initialise temporary dictionary
      call tempDict % init(4,1)

      ! Build approperiate surface
      if (self % r(i) == ZERO ) then
        call tempDict % store('type','infSurf')
        call tempDict % store('id', 1 )

      else
        call tempDict % store('type','zCylinder')
        call tempDict % store('id',1)
        call tempDict % store('origin',[ZERO, ZERO, ZERO])
        call tempDict % store('radius', self % r(i))

      end if

      ! Put surface on the shelf
      allocate(tempSurf, source = new_surface(tempDict))
      call sShelf % getOrAdd(tempSurf, surfIDs(i), self % data(i) % surfIdx)
      call tempDict % kill()
    end do

    ! Build innermost cell
    call tempCell % init( [-surfIDs(1)], 1, sShelf)
    call cShelf % getOrAdd(tempCell, cellId, self % data(1) % cellIdx)

    ! Build rest of the cells
    do i=2,size(self % r)
      call tempCell % init([surfIDs(i-1), -surfIDs(i)], 1, sShelf)
      call cShelf % getOrAdd(tempCell, cellId, self % data(i) % cellIdx)

    end do
  end subroutine buildAnnuli

  !!
  !! Translate material names to fillVector
  !!
  subroutine makeFillVector(self, fillVEctor, matNames ,materials)
    class(pinUniverse), intent(inout)                      :: self
    integer(shortInt),dimension(:),allocatable,intent(out) :: fillVector
    character(nameLen),dimension(:),intent(in)             :: matNames
    class(nuclearData), intent(in)                         :: materials
    integer(shortInt)                                      :: i, fill, L
    character(nameLen)                                     :: tempName
    character(100), parameter :: Here ='makeFillVector (pinUniverse_class.f90)'

    ! Double check that size of radii match size of matNames
    if( size(self % r) /= size(matNames)) call fatalError(Here,'matNames not matching self % r')

    ! Allocate fill vector
    allocate( fillVector( size(matNames)))

    ! Loop over all matNames
    do i=1, size(matNames)
      ! Copy material name, adjust left and find length of name
      tempName = adjustl(matNames(i))
      L = len_trim(tempName)

      ! Verify if name is a universe idd
      if ( tempName(1:2) == 'u<' .and. tempName(L:L) == '>') then
        ! Read universe id into fill
        tempName = tempName(3:L-1)
        read(tempName,*) fill
        fill = -fill

      else
        fill = materials % getIdx(tempName)

      end if

      ! Load fill into vector
      fillVector(i) = fill

    end do
  end subroutine makeFillVector


end module pinUniverse_class
