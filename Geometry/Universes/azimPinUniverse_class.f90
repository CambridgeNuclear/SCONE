module azimPinUniverse_class

  use numPrecision
  use universalVariables, only : INF, targetNotFound
  use genericProcedures,  only : fatalError, numToChar, swap
  use dictionary_class,   only : dictionary
  use coord_class,        only : coord
  use charMap_class,      only : charMap
  use surfaceShelf_class, only : surfaceShelf
  use cylinder_class,     only : cylinder
  use plane_class,        only : plane
  use cell_inter,         only : cell
  use cellShelf_class,    only : cellShelf
  use universe_inter,     only : universe, kill_super => kill, charToFill
  implicit none
  private

  ! Parameters
  ! Are public for use in unit tests
  integer(shortInt), parameter, public :: MOVING_IN = -1, MOVING_OUT = -2, &
                                          MOVING_ANTI = -3, MOVING_CLOCK = -4, &
                                          MOVING_CLOCK_BACK = -5, &
                                          MOVING_CLOCK_FORWARD = -6

  !!
  !! Universe that represents a single pin
  !! A version of the pinUniverse which is equi-azimuthally divided.
  !!
  !! Is composed from co-centring cylinders divided azimuthally by planes. 
  !! In the most central cylinder and first azimuthal segment, the cell will
  !! have an ID of 1. This increases by 1 while proceeding across azimuthal
  !! segments. Proceeding radially outwards, the ID is incremented by the 
  !! number of azimuthal regions.
  !!
  !! This has been generalised to allow different numbers of azimuthal divisions
  !! in each radial division, e.g., a pin cell with 4 azimuthal divisions in the
  !! fuel but 8 azimuthal divisions in the moderator
  !!
  !! Simplified input for uniform azimuthal division
  !! Sample Dictionary Input:
  !!   azimPinUni {
  !!     id 7;
  !!     type azimPinUniverse;
  !!     naz 4;
  !!     #origin (1.0 0.0 0.1);#
  !!     #rotation (30.0 0.0 0.0);#
  !!     radii (3.0 4.5 0.0 1.0 );
  !!     fills (u<3> void clad u<4>);
  !!   }
  !!
  !! Input for non-uniform azimuthal division
  !! Sample Dictionary Input:
  !!   azimPinUni {
  !!     id 7;
  !!     type azimPinUniverse;
  !!     nazR (4 8);
  !!     #origin (1.0 0.0 0.1);#
  !!     #rotation (30.0 0.0 0.0);#
  !!     radii (3.0 0.0 );
  !!     fills (u<3> water );
  !!   }
  !!
  !!  naz corresponds to the number of azimuthal regions produced. Must be a multiple of 2.
  !!  Takes origin at 0 degrees, i.e., the centre of the first azimuthal slice.
  !!  There must be 0.0 entry, which indicates outermost annulus (infinite radius).
  !!
  !!  Alternatively, nazR is the number of azimuthal divisions per radial region, from
  !!  the centre radiating outwards.
  !!
  !!  `fills` and `radii` are given as pairs by position in the input arrays. Thus, fills
  !!  are sorted together with the `radii`. As a result, in the example, local cells 1 to 4
  !!  are filled with u<4>, cell 5 to 8 with u<3> etc.
  !!
  !!  !!!!!
  !!  TODO: Just for the moment, there are no azimuthally different fills. This should be remedied!
  !!  !!!!!
  !!
  !! Public Members:
  !!  nAz     -> Number of azimuthal regions in each radial ring
  !!  r_sqr   -> Array of radius^2 for each annulus
  !!  theta   -> Array of azimuthal boundary angles in radians. 
  !!  annuli  -> Array of cylinder surfaces that represent diffrent annuli
  !!  planes  -> Array of planes for providing azimuthal division
  !!  normals -> Array of plane normals for convenience
  !!
  !! Interface:
  !!   universe interface
  !!
  type, public, extends(universe) :: azimPinUniverse
    private
    integer(shortInt), dimension(:), allocatable :: nAz
    real(defReal), dimension(:), allocatable     :: r_sq
    real(defReal), dimension(:), allocatable     :: theta
    type(cylinder), dimension(:), allocatable    :: annuli
    type(plane), dimension(:), allocatable       :: planes
    real(defReal), dimension(:,:), allocatable   :: normals
  contains
    ! Superclass procedures
    procedure :: init
    procedure :: kill
    procedure :: findCell
    procedure :: distance
    procedure :: cross
    procedure :: cellOffset
  end type azimPinUniverse

contains

  !!
  !! Initialise Universe
  !!
  !! See universe_inter for details.
  !!
  !! Errors:
  !!   fatalError for invalid input
  !!
  subroutine init(self, fill, dict, cells, surfs, mats)
    class(azimPinUniverse), intent(inout)                     :: self
    integer(shortInt), dimension(:), allocatable, intent(out) :: fill
    class(dictionary), intent(in)                             :: dict
    type(cellShelf), intent(inout)                            :: cells
    type(surfaceShelf), intent(inout)                         :: surfs
    type(charMap), intent(in)                                 :: mats
    integer(shortInt)                             :: id, idx, N, i, j, r, nAzTotal
    real(defReal), dimension(:), allocatable      :: radii, temp
    integer(shortInt), dimension(:), allocatable  :: tempInt
    character(nameLen), dimension(:), allocatable :: fillNames
    real(defReal)                                 :: dTheta, theta0
    character(100), parameter :: Here = 'init (azimPinUniverse_class.f90)'

    ! Load basic data
    call dict % get(id, 'id')
    if (id <= 0) call fatalError(Here, 'Universe ID must be +ve. Is: '//numToChar(id))
    call self % setId(id)
    
    ! Load radii and fill data
    call dict % get(radii, 'radii')
    call dict % get(fillNames, 'fills')

    ! Check values
    if (size(radii) /= size(fillNames)) then
      call fatalError(Here, 'Size of radii and fills does not match')

    else if (any(radii < ZERO)) then
      call fatalError(Here, 'Found -ve value of radius.')

    end if

    ! Load azimuthal division
    if (dict % isPresent('naz') .and. dict % isPresent('nazR')) then
      call fatalError(Here,'Cannot have both a naz and nazR entry')

    ! Only one azimuthal discretisation
    elseif (dict % isPresent('naz')) then
      allocate(self % nAz(size(radii)))
      call dict % get(N, 'naz')
      self % nAz = N

    ! Variable azimuthal discretisation
    elseif (dict % isPresent('nazR')) then
      call dict % get(tempInt, 'nazR')
      if (size(tempInt) /= size(radii)) call fatalError(Here,'Number of radial regions is not consistent '//&
              'between nazR and radii')
      allocate(self % nAz(size(tempInt)))
      self % nAz = tempInt
    else
      call fatalError(Here,'Must have either a naz or nazR entry')
    end if
    if (any(self % nAz < 2)) call fatalError(Here,'Number of azimuthal regions must be 2 or more')  

    ! Use binary logic to check if nAz is a power of 2
    do r = 1, size(self % naz)
      if (IAND(self % nAz(r), self % nAz(r) - 1) /= 0) call fatalError(Here, 'Number of azimuthal regions must be a power of 2')
    end do

    ! Load origin
    if (dict % isPresent('origin')) then
      call dict % get(temp, 'origin')

      if (size(temp) /= 3) then
        call fatalError(Here, 'Origin must have size 3. Has: '//numToChar(size(temp)))
      end if
      call self % setTransform(origin=temp)

    end if

    ! Load rotation
    if (dict % isPresent('rotation')) then
      call dict % get(temp, 'rotation')

      if (size(temp) /= 3) then
        call fatalError(Here, '3 rotation angles must be given. Has: '//numToChar(size(temp)))
      end if
      call self % setTransform(rotation=temp)
    end if


    ! Sort radii with selection sort
    ! Start with value 0.0 that represents outermost element
    ! Change 0.0 to infinity
    N = size(radii)
    idx = minloc(radii, 1)
    if (radii(idx) /= ZERO) call fatalError(Here, 'Did not find outermost element with radius 0.0.')
    call swap( radii(idx), radii(N))
    call swap( fillNames(idx), fillNames(N))
    radii(N) = INF * 1.1_defReal

    do i = N-1,1,-1
      idx = maxloc(radii(1:i), 1)
      call swap( radii(idx), radii(i))
      call swap( fillNames(idx), fillNames(i))
    end do

    ! Check for duplicate values of radii
    do i = 1, N-1
      if (radii(i) == radii(i+1)) then
        call fatalError(Here, 'Duplicate value of radius: '//numToChar(radii(i)))
      end if
    end do

    ! Load data & Build cylinders
    self % r_sq = radii * radii

    allocate(self % annuli(N))
    do i = 1, N
      call self % annuli(i) % build(id=1, origin=[ZERO, ZERO, ZERO], &
                                    type='zCylinder', radius=radii(i) )
    end do

    ! Load data and build planes
    allocate(self % theta(sum(self % nAz)))
    allocate(self % planes(sum(self % nAz)/2))
    allocate(self % normals(sum(self % nAz)/2,2))

    nAzTotal = 0
    do r = 1, N
    
      dTheta = TWO_PI / self % nAz(r)
      theta0 = HALF * dTheta
      do i = 1, self % nAz(r)
        self % theta(nAzTotal + i) = theta0 
        theta0 = theta0 + dTheta
      end do

      ! Build the planes, rotated at angle theta from the line x = 0
      do i = 1, self % nAz(r)/2
        self % normals(nAzTotal/2 + i,1) = -sin(self % theta(nAzTotal + i)) 
        self % normals(nAzTotal/2 + i,2) = cos(self % theta(nAzTotal + i))
        call self % planes(nAzTotal/2 + i) % build(id=1, &
                norm = [self % normals(nAzTotal/2 + i,1), self % normals(nAzTotal/2 + i,2), ZERO], offset=ZERO)
      end do

      nAzTotal = nAzTotal + self % nAz(r)
    end do

    ! Create fill array
    allocate(fill(nAzTotal))

    nAzTotal = 0
    do i = 1, N
      do j = 1, self % nAz(i)
        fill(nAzTotal + j) = charToFill(fillNames(i), mats, Here)
      end do
      nAzTotal = nAzTotal + self % nAz(i)
    end do

  end subroutine init

  !!
  !! Find local cell ID given a point
  !!
  !! See universe_inter for details.
  !!
  subroutine findCell(self, localID, cellIdx, r, u)
    class(azimPinUniverse), intent(inout)   :: self
    integer(shortInt), intent(out)          :: localID
    integer(shortInt), intent(out)          :: cellIdx
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    real(defReal)                           :: r_sq, theta, mul, planeDir
    integer(shortInt)                       :: aIdx, rIdx, pIdx, baseAIdx

    r_sq = r(1)*r(1) + r(2)*r(2)
    theta = atan2(r(2),r(1))
    ! Correct theta for negative values
    if (theta < ZERO) theta = TWO_PI + theta

    cellIdx = 0

    ! Need to include surface tolerance. Determine multiplier by direction
    if ( r(1)*u(1) + r(2)*u(2) >= ZERO) then
      mul = -ONE
    else
      mul = ONE
    end if

    ! Find local cell
    ! Start by finding annular region
    do rIdx = 1, size(self % r_sq)
      if( r_sq < self % r_sq(rIdx) + mul * self % annuli(rIdx) % surfTol() ) exit
    end do
    ! If reached here without exiting, rIdx = size(self % r_sq) + 1

    ! Find base azimuthal index given radial zone
    baseAIdx = 0
    if (rIdx > 1) baseAIdx = sum(self % nAz(1:(rIdx-1)))

    ! Find azimuthal segment
    do aIdx = 1, self % nAz(rIdx)
      if (aIdx > self % nAz(rIdx)/2) then
        pIdx = aIdx - self % nAz(rIdx)/2
        planeDir = -ONE
      else
        pIdx = aIdx
        planeDir = ONE
      end if
      ! Surface tolerance multiplier determined by relative direction of particle
      ! and theta
      if (planeDir*self % normals(baseAIdx/2 + pIdx,1)*u(1) + planeDir*self % normals(baseAIdx/2 + pIdx,2)*u(2) >= ZERO) then
        mul = -ONE
      else
        mul = ONE
      end if
      if (theta < self % theta(baseAIdx + aIdx) + mul * self % planes(baseAIdx/2 + pIdx) % surfTol() ) exit
    end do
    ! If exceeded the search, theta is <2pi but greater than the largest theta.
    ! Therefore, it lies in the negative theta portion of the first segment.
    if (aIdx == self % nAz(rIdx) + 1) aIdx = 1
    localID = aIdx + baseAIdx

  end subroutine findCell

  !!
  !! Return distance to the next boundary between local cells in the universe
  !!
  !! See universe_inter for details.
  !!
  !! Errors:
  !!   fatalError is localID is invalid
  !!
  subroutine distance(self, d, surfIdx, coords)
    class(azimPinUniverse), intent(inout) :: self
    real(defReal), intent(out)            :: d
    integer(shortInt), intent(out)        :: surfIdx
    type(coord), intent(in)               :: coords
    real(defReal)                         :: d_out_annul, d_in_annul
    real(defReal)                         :: d_plus, d_minus, dAz
    integer(shortInt)                     :: id, aIdx, rIdx, i, searchIdxP, searchIdxM
    integer(shortInt)                     :: sense, minus_sense, plus_sense, baseIdx
    character(100), parameter :: Here = 'distance (azimPinUniverse_class.f90)'

    ! Get local id
    id = coords % localID

    if (id < 1 .or. id > sum(self % nAz)) then
      call fatalError(Here, 'Invalid local ID: '//numToChar(id))
    end if

    ! Identify annulus index and azimuthal index
    do i = 1, size(self % annuli)
      id = id - self % nAz(i)
      if (id < 1) then
        rIdx = i
        aIdx = id + self % nAz(i) 
        exit
      end if
    end do

    baseIdx = 0
    if (rIdx > 1) baseIdx = sum(self % nAz(1:(rIdx-1)))/2

    ! Check distance to annuli
    ! Outer distance
    if (rIdx > size(self % r_sq)) then
      d_out_annul = INF
    else
      d_out_annul = self % annuli(rIdx) % distance(coords % r, coords % dir)
    end if

    ! Inner distance
    if (rIdx == 1) then
      d_in_annul = INF
    else
      d_in_annul = self % annuli(rIdx-1) % distance(coords % r, coords % dir)
    end if

    ! Select distance and surface
    if ( d_in_annul < d_out_annul) then
      surfIdx = MOVING_IN
      d = d_in_annul

    else
      surfIdx = MOVING_OUT
      d = d_out_annul
    end if

    ! Check distance to azimuthal planes

    ! Check whether in the first or second half of azimuthal segments
    ! If in second half, find the same azimuthal indices as in the first half
    searchIdxP = aIdx + baseIdx
    if (aIdx > self % nAz(rIdx)/2) searchIdxP = searchIdxP - self % nAz(rIdx)/2

    ! Set default senses for which cell will be entered
    minus_sense = MOVING_CLOCK
    plus_sense  = MOVING_ANTI

    ! Check for first or last cells circling round
    if (aIdx == 1) then
      minus_sense = MOVING_CLOCK_BACK
    else if (aIdx == self % nAz(rIdx)) then
      plus_sense = MOVING_CLOCK_FORWARD
    end if

    ! If in the first cell or nAz/2 + 1 cell, find correct plane
    if (searchIdxP == baseIdx + 1) then
      searchIdxM = self % nAz(rIdx)/2 + baseIdx
    else
      searchIdxM = searchIdxP - 1
    end if

    ! Identify which two planes (or only one if nAz = 2)
    ! Check to see if in second half of azimuthal segments
    if (self % nAz(rIdx) > 2) then
      d_plus  = self % planes(searchIdxP) % distance(coords % r, coords % dir)
      d_minus = self % planes(searchIdxM) % distance(coords % r, coords % dir)
    else
      d_plus  = self % planes(baseIdx + 1) % distance(coords % r, coords % dir)
      d_minus = INF
    end if

    ! Choose minimum azimuthal
    if (d_plus < d_minus) then
      dAz = d_plus
      sense = plus_sense
    else
      dAz = d_minus
      sense = minus_sense
    end if
    
    ! Compare radial distance with azimuthal
    if (dAz < d) then
      surfIdx = sense
      d = dAz
    end if

  end subroutine distance

  !!
  !! Cross between local cells
  !!
  !! See universe_inter for details.
  !!
  !! Errors:
  !!   fatalError if surface from distance is not MOVING_IN or MOVING_OUT
  !!
  subroutine cross(self, coords, surfIdx)
    class(azimPinUniverse), intent(inout) :: self
    type(coord), intent(inout)            :: coords
    integer(shortInt), intent(in)         :: surfIdx
    integer(shortInt)                     :: aIdx, i
    character(100), parameter :: Here = 'cross (azimPinUniverse_class.f90)'
    
    ! Need radial region to work out how much to increment the clock by
    ! Identify annulus index and azimuthal index
    aIdx = coords % localID
    do i = 1, size(self % annuli)
      aIdx = aIdx - self % nAz(i)
      if (aIdx < 1) exit
    end do

    ! Need to determine whether completing a circle, i.e., moving from
    ! the last azimuthal segment to the first and vice versa
    if (surfIdx == MOVING_CLOCK) then
      coords % localID = coords % localID - 1

    else if (surfIdx == MOVING_ANTI) then
      coords % localID = coords % localID + 1

    else if (surfIdx == MOVING_CLOCK_BACK) then
      coords % localID = coords % localID + self % nAz(i) - 1

    else if (surfIdx == MOVING_CLOCK_FORWARD) then
      coords % localID = coords % localID - self % nAz(i) + 1

    ! Need to be replaced with find subroutines in the general case, sadly
    else if (surfIdx == MOVING_IN) then
      if (self % nAz(i) == self % nAz(i-1)) then
        coords % localID = coords % localID - self % nAz(i)
      else
        call self % findCell(coords % localID, coords % cellIdx, coords % r, coords % dir)
      end if
    else if (surfIdx == MOVING_OUT) then
      if (self % nAz(i) == self % nAz(i+1)) then
        coords % localID = coords % localID + self % nAz(i)
      else
        call self % findCell(coords % localID, coords % cellIdx, coords % r, coords % dir)
      end if
    else
      call fatalError(Here, 'Unknown surface memento: '//numToChar(surfIdx))

    end if

  end subroutine cross

  !!
  !! Return offset for the current cell
  !!
  !! See universe_inter for details.
  !!
  function cellOffset(self, coords) result (offset)
    class(azimPinUniverse), intent(in) :: self
    type(coord), intent(in)            :: coords
    real(defReal), dimension(3)        :: offset

    ! There is no cell offset
    offset = ZERO

  end function cellOffset

  !!
  !! Return to uninitialised state
  !!
  subroutine kill(self)
    class(azimPinUniverse), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Kill local
    if(allocated(self % r_sq)) deallocate(self % r_sq)
    if(allocated(self % annuli)) deallocate(self % annuli)
    if(allocated(self % theta)) deallocate(self % theta)
    if(allocated(self % planes)) deallocate(self % planes)
    if(allocated(self % normals)) deallocate(self % normals)
    if(allocated(self % nAz)) deallocate(self % nAz)

  end subroutine kill

end module azimPinUniverse_class
