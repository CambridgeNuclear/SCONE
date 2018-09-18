module latUniverse_class

  use numPrecision
  use universalVariables
  use genericProcedures,   only : fatalError, targetNotFound, numToChar
  use vector_class,        only : vector
  use dictionary_class,    only : dictionary
  use coord_class,         only : coord
  use maps_class,          only : intMap
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
  interface latUniverse
    module procedure latUniverse_fromDict
  end interface

  !!
  !! Regular rectangular lattice in 2D or 3D
  !!
  !! Enumeration is from lower, bottom left corner inside the representation of the latticle
  !! z = 1  z = 2   z = 3
  !!  3  4   7  8   11 12
  !!  1  2   5  6   9  10
  !!
  type, public, extends(universe) :: latUniverse
    private
    integer(shortInt)               :: maxDim = 0
    real(defReal), dimension(3)     :: pitch  = 0.0
    integer(shortInt), dimension(3) :: sizeN  = 0
    real(defReal), dimension(3)     :: corner
    real(defReal), dimension(3)     :: half_bar

    integer(shortInt)               :: pCellIdx   ! Pitch cell index
    integer(shortInt)               :: oCellIdx   ! Out cell index
    integer(shortInt)               :: outLocalID

  contains
    ! Build procedures
    procedure :: init
    procedure :: cellIdx

    ! Runtime procedures
    procedure :: findCell
    procedure :: distance
    procedure :: cross
    procedure :: cellOffset

  end type latUniverse

contains
  !!
  !! Initialise pin Universe
  !! Creates only geometry definition. fillVector is build in constructor function
  !!
  subroutine init(self, offset, id, pitch, sizeN, sShelf, cShelf)
    class(latUniverse), intent(inout)           :: self
    real(defReal), dimension(:), intent(in)     :: offset
    integer(shortInt), intent(in)               :: id
    real(defReal), dimension(:), intent(in)     :: pitch
    integer(shortInt), dimension(:), intent(in) :: sizeN
    type(surfaceShelf), intent(inout)           :: sShelf
    type(cellShelf),intent(inout)               :: cShelf
    integer(shortInt)                           :: maxDim
    type(dictionary)                            :: cellDict, outDict
    class(surface), allocatable                 :: cellSurf, outSurf
    type(cell)                                  :: pitchCell, outCell
    real(defReal), dimension(3)                 :: cellHW, outHW, offsetTemp
    character(nameLen)                          :: surfType
    integer(shortInt)                           :: cellSId, outSId, dummy
    character(100), parameter :: Here = 'init (latUniverse_class.f90)'

    ! Verify size and choose whether lattice is 2D or 3D
    maxDim = size(offset)
    if( size(pitch) /= maxDim) call fatalError(Here,"Pitch doesn't have size : "//numToChar(maxDim))
    if( size(sizeN) /= maxDim) call fatalError(Here,"SizeN doesn't have size : "//numToChar(maxDim))

    ! Check max dim
    if (2 > maxDim .or. maxDim > 3) call fatalError(Here,"Dimension should be 2 or 3")

    ! Check invalid pitch
    if (any(pitch < surface_tol)) then
      call fatalError(Here, "Pitch must be greater then surface tolerance")
    end if

    ! Check invalid size
    if (any(sizeN < 1)) then
      call fatalError(Here, "Lattice size in every dimension must be > then 1")
    end if

    ! Load id and offset
    call self % setId(id)
    if (maxDim == 2) then
      offsetTemp = [offset, ZERO]
    else
      offsetTemp = offset
    end if
    call self % setOffset(offsetTemp)

    ! Load lattice parameters
    self % maxDim          = maxDim
    self % pitch(1:maxDim) = pitch
    self % sizeN(1:maxDim) = sizeN

    ! Select surface type and halfwidths based on dimension
    select case(maxDim)
      case(2)
        surfType = 'zSquareCylinder'
        cellHW  =  [pitch, ZERO]/2.0
        outHW   =  [pitch, ZERO] * [sizeN, 0] / 2.0

      case(3)
        surfType = 'box'
        cellHW  =  pitch / 2.0
        outHW   =  pitch * sizeN / 2.0

      case default
        call fatalError(Here,'This should never happen. God must have been decleard integer!')

    end select

    ! Define cell surface
    call cellDict % init(4,1)
    call cellDict % store('type',surfType)
    call cellDict % store('id',1)
    call cellDict % store('origin',[ZERO, ZERO, ZERO])
    call cellDict % store('halfwidth',cellHW)

    ! Define outside surface
    call outDict % init(4,1)
    call outDict % store('type',surfType)
    call outDict % store('id',1)
    call outDict % store('origin',offsetTemp)
    call outDict % store('halfwidth',outHW)

    ! Add surfaces to the shelf
    allocate( cellSurf, source = new_surface(cellDict))
    allocate( outSurf,  source = new_surface(outDict))

    call sShelf % getOrAdd(cellSurf, cellSId, dummy)
    call sShelf % getOrAdd(outSurf,  outSId,  dummy)

    ! Create lattice cell and outside cell & add to shelf
    call pitchCell % init([-cellSId],1,sShelf)
    call outCell   % init([ outSId], 1,sShelf)

    call cShelf % getOrAdd(pitchCell, dummy, self % pCellIdx)
    call cShelf % getOrAdd(outCell,   dummy, self % oCellIdx)

    ! Save local id of outside cell
    self % outLocalID = product(sizeN, sizeN > 0) + 1

    ! Save position of the bottom(Z) bottom(Y) left corner
    self % corner = - outHW

    ! Save surface tolerance corrested cell halfwidth
    self % half_bar = self % pitch * HALF - surface_tol

  end subroutine init

  !!
  !! Returns an initialised instance of a lattice universe
  !! Returns fillVector as well with +ve entries being materialIDXs and -ve fill universes IDs
  !! Provisionally provide nuclearData *** REPLACE WITH CHAR-INT MAP LATER FOR DECOUPLING
  !! *** TODO: IMPLEMENT MAP REPEAT FEATURE
  !!
  function latUniverse_fromDict(fillVector, dict, cShelf, sShelf, cellFillMap, materials) result (new)
    integer(shortInt),dimension(:),allocatable,intent(out) :: fillVector
    class(dictionary), intent(in)                          :: dict
    type(cellShelf), intent(inout)                         :: cShelf
    type(surfaceShelf), intent(inout)                      :: sShelf
    type(intMap), intent(in)                               :: cellFillMap
    class(nuclearData), intent(in)                         :: materials
    type(latUniverse)                                      :: new
    integer(shortInt)                                      :: id, maxDim
    integer(shortInt)                                      :: mapSize
    real(defReal), dimension(:), allocatable               :: offset, pitch
    character(nameLen)                                     :: rep, padMat
    integer(shortInt), dimension(:), allocatable           :: sizeN
    character(100), parameter :: Here = 'latUniverse_fromDict (latUniverse_class.f90)'

    ! Read universe data
    call dict % get(id,'id')
    call dict % get(offset,'origin')
    call dict % get(pitch,'pitch')
    call dict % get(sizeN,'shape')
    call dict % get(fillVector,'map')
    call dict % get(padMat,'padMat')

    call dict % getOrDefault(rep,'repeat','none')

    ! Check maximum dimension
    maxDim = size(sizeN)
    if (maxDim == 3) then
      if (sizeN(3) == 0 ) maxDim = 2
    end if

    ! Initialise universe
    call new % init( offset(1:maxDim), id, pitch(1:maxDim), sizeN(1:maxDim), sShelf, cShelf)

    ! Check fill vector
    if (any(fillVector < 1)) call fatalError(Here,'Invalid universe IDs in lattice map (uniId < 1)')

    ! Check size of lattice map. mapSize -> required size of universe map
    mapSize = product(sizeN(1:maxDim))
    if (size(fillVector) /= mapSize) then
      call fatalError(Here,"Universe map has size of: "//numToChar(size(fillVector)) &
                           //" should have size of: " // numToChar(mapSize))
    end if

    ! Build fill vector
    fillVector = -fillVector
    fillVector = [fillVector, materials % getIdx(padMat)]

  end function latUniverse_fromDict

  !!
  !! Return cell index given localID of a cell
  !!
  function cellIdx(self,localId)
    class(latUniverse), intent(in) :: self
    integer(shortInt), intent(in)  :: localID
    integer(shortInt)              :: cellIdx
    character(100), parameter :: Here = 'cellIdx (latUniverse_class.f90)'

    if (localID < 1 .or. localID > self % outLocalID) then
      call fatalError(Here,'Provided local ID ' // numToChar(localID) // ' is Not present')
    end if

    if (localID == self % outLocalID) then
      cellIdx = self % oCellIdx

    else
      cellIdx = self % pCellIdx

    end if
  end function cellIdx

  !!
  !! Using the coordinates it finds a localID & cellIDx inside the universe
  !!
  subroutine findCell(self, coords, cShelf, sShelf)
    class(latUniverse), intent(in)  :: self
    type(coord),intent(inout)       :: coords
    type(cellShelf), intent(in)     :: cShelf
    type(surfaceShelf), intent(in)  :: sShelf
    integer(shortInt)               :: maxDim, i, inc
    real(defReal),dimension(3)      :: r_bar,u
    logical(defBool)                :: outside, atSurface, goingOut
    integer(shortInt), dimension(3) :: ijk

    maxDim = self % maxDim
    ijk = 1

    r_bar = coords % r - self % corner
    u     = coords % dir

    ! Calculate position inside the lattice
    ijk(1:maxDim) = floor( r_bar(1:maxDim) / self % pitch(1:maxDim)) + 1

    ! Calculate coordinates wrt middle of the current cell
    r_bar(1:maxDim) = r_bar(1:maxDim) - (ijk(1:maxDim) - HALF) * self % pitch(1:maxDim)

    ! Check if position is within surface tolerance and if it is moving outside
    ! Use the fact that it is moving outside if u(i) and r(i) have the same sign
    do i=1,maxDim
      atSurface = abs(r_bar(i)) > self % half_bar(i)
      goingOut  = u(i) * r_bar(i) > ZERO

      if (atSurface .and. goingOut ) then
        ! Set increment to 1 if u is +ve and -1 if u is -ve
        inc = 1
        if ( u(i) < 0 ) inc = -1

        ! Update ijk index and local coordinates
        ijk(i) = ijk(i) + inc
        r_bar(i) = r_bar(i) - self % pitch(i) * inc
      end if
    end do

    ! Check if point is outside
    outside = .false.
    do i=1,maxDim
      outside = outside .or. (ijk(i) <= 0 .or. ijk(i) > self % sizeN(i))

    end do

    ! Set localID and cellIdx
    if (outside) then
      coords % localID = self % outLocalID
      coords % cellIdx = self % oCellIdx

    else
      coords % localID = ijk(1) + self % sizeN(1) * (ijk(2)-1 + self % sizeN(2) * (ijk(3)-1) )
      coords % cellIdx = self % pCellIdx

      ! Shift local coords to correspond to centre of a pitch cell
      coords % r(1:maxDim) = r_bar(1:maxDim)

    end if

  end subroutine findCell

  !!
  !! Returns distance to the next cell boundary in the universe
  !! Returns surfIdx of the surface beeing X-ed
  !! surfIdx can be < 0 (invalid)
  !! Invalid surfIdx indicate a surface that is private to the universe(not present on sShelf)
  !! Example of a private surface is face of a lattice cell
  !!
  subroutine distance(self, dist, surfIdx, coords ,cShelf, sShelf)
    class(latUniverse), intent(in) :: self
    real(defReal), intent(out)     :: dist
    integer(shortInt), intent(out) :: surfIdx
    type(coord), intent(in)        :: coords
    type(cellShelf), intent(in)    :: cShelf
    type(surfaceShelf),intent(in)  :: sShelf
    real(defReal),dimension(3)     :: bound
    real(defReal),dimension(3)     :: r, u
    type(vector)                   :: rv, uv
    real(defReal)                  :: testDist
    integer(shortInt)              :: maxDim, i, xDir

    if(coords % localID /= self % outLocalID) then

      ! Load maximum dimension & position & direction
      maxDim = self % maxDim
      r = coords % r
      u = coords % dir

      ! Find bounds of the pitch cell
      bound = sign(self % pitch * HALF, u)

      ! Find minimum distance and direction
      dist = INFINITY

      do i=1,maxDim
        ! Calculate distance to the boundary
        if ( u(i) /= ZERO) then
          testDist = (bound(i) - r(i)) / u(i)

        else
          testDist = INFINITY

        end if

        testDist = max(ZERO,testDist)

        ! Save shortest distance
        if (testDist < dist) then
          dist = testDist
          xDir = i
        end if

      end do

      ! Find code of a face beeing X-ed
      surfIdx = xDir * 2 - 1
      if ( u(xDir) < ZERO ) surfIdx = surfIdx +1
      surfIdx = -surfIdx
    else
      ! Get position and direction
      rv = coords % r
      uv = coords % dir

      call cShelf % shelf(self % oCellIdx) % distance(dist, surfIdx, rv, uv, sShelf)
    end if
  end subroutine distance

  !!
  !! Perform crossing inside the universe from current cell to next cell
  !! Assumes coords are at the surface being crossed (within surface_tol)
  !! IT DOES NOT PERFORM CHECK IF THE ABOVE ASSUMPTION IS VALID! (Be careful - MAK)
  !!
  subroutine cross(self, coords, surfIdx, cShelf, sShelf)
    class(latUniverse), intent(in)  :: self
    type(coord),intent(inout)       :: coords
    integer(shortInt), intent(in)   :: surfIdx
    type(cellShelf), intent(in)     :: cShelf
    type(surfaceShelf), intent(in)  :: sShelf
    integer(shortInt), dimension(3) :: ijk
    integer(shortInt)               :: localId, inc, dir, maxDim

    ! Set local cell to outside. Save localID outside bounds to reconstruct ijk
    localID = coords % localID

    if (coords % localID /= self % outLocalID) then
      ! Find direction beeing X-ed (Note that surfIdx is -ve)
      dir = (1-surfIdx) / 2

      ! Set increment to 1 if u is +ve and -1 if u is -ve
      inc = 1
      if ( coords % dir(dir) < 0 ) inc = -1

      ! Increment location
      coords % r(dir) = coords % r(dir) - self % pitch(dir) * inc
      coords % localID = coords % localID + inc * product(self % sizeN(1:dir-1))

      ! If partcile has leaked change to frame wrt centre of the lattice
      if (coords % localID < 1 .or. coords % localID >= self % outLocalID) then
        ! Set localID to outside cell
        coords % localID = self % outLocalID

        ! Load maximum dimenstion
        maxDim = self % maxDim

        ! Return to co-ordinates aligned with centre of lattice
        ijk(1) = modulo(localID-1, self % sizeN(1)) + 1
        localID = (localID - ijk(1)) / self % sizeN(1)
        ijk(2) = modulo(localID, self % sizeN(2)) + 1
        localID = (localID - ijk(2) +1) / self % sizeN(2)
        ijk(3) = localID + 1

        ! Reapply crossing on ijk
        ijk(dir) = ijk(dir) + inc

        coords % r(1:maxDim) = self % corner(1:maxDim) + &
                             (ijk(1:maxDim) - HALF) * self % pitch(1:maxDim) + coords % r(1:maxDim)
      end if

    else
      call self % findCell(coords,cShelf,sShelf)

    end if
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
    class(latUniverse), intent(in) :: self
    type(coord),intent(inout)      :: coords
    real(defReal),dimension(3)     :: offset

    offset = [ZERO, ZERO, ZERO]

  end function cellOffset
    
end module latUniverse_class
