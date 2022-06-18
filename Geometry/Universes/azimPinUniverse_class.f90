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
  !! Is composed from co-centring cylinders. Central cell has local ID 1 and the ID
  !! increases with subsequent rings.
  !!
  !! Sample Dictionary Input:
  !!   pinUni {
  !!     id 7;
  !!     type azimPinUniverse;
  !!     #naz 4#
  !!     #origin (1.0 0.0 0.1);#
  !!     #rotation (30.0 0.0 0.0);#
  !!     radii (3.0 4.5 0.0 1.0 );
  !!     fills (u<3> void clad u<4>);
  !!   }
  !!
  !!  naz corresponds to the number of azimuthal regions produced. Must be a multiple
  !!  of 2.
  !!  Takes origin at 0 degrees, i.e., the centre of the first azimuthal slices.
  !!  There must be 0.0 entry, which indicates outermost annulus (infinite radius).
  !!  `fills` and `radii` are given as pairs by position in the input arrays. Thus, fills
  !!  are sorted together with the `radii`. As a result, in the example, local cell 1 is
  !!  filled with u<4>, cell 2 with u<3> etc.
  !!
  !!  !!!!!
  !!  Just for the moment, there are no azimuthally different fills. This will be remedied!!
  !!  !!!!!
  !!
  !! Public Members:
  !!  nAz    -> Number of azimuthal regions
  !!  r_sqr  -> Array of radius^2 for each annulus
  !!  theta  -> Array of azimuthal boundary angles in radians. 
  !!  annuli -> Array of cylinder surfaces that represent diffrent annuli
  !!
  !! Interface:
  !!   universe interface
  !!
  type, public, extends(universe) :: azimPinUniverse
    private
    integer(shortInt)                          :: nAz
    real(defReal), dimension(:), allocatable   :: r_sq
    real(defReal), dimension(:), allocatable   :: theta
    type(cylinder), dimension(:), allocatable  :: annuli
    type(plane), dimension(:), allocatable     :: planes
    real(defReal), dimension(:,:), allocatable :: normals
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
    integer(shortInt)                             :: id, idx, N, i, j
    real(defReal), dimension(:), allocatable      :: radii, temp
    character(nameLen), dimension(:), allocatable :: fillNames
    real(defReal)                                 :: dTheta, theta0
    character(100), parameter :: Here = 'init (azimPinUniverse_class.f90)'

    ! Load basic data
    call dict % get(id, 'id')
    if (id <= 0) call fatalError(Here, 'Universe ID must be +ve. Is: '//numToChar(id))
    call self % setId(id)

    ! Load azimuthal division
    call dict % get(self % naz, 'naz')
    if (self % nAz < 2) call fatalError(Here,'Number of azimuthal regions must be 2 or more')

    ! Use binary logic to check if nAz is a power of 2
    if (IAND(self % nAz, self % nAz - 1) /= 0) call fatalError(Here, 'Number of azimuthal regions must be a multiple of 2')

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
        call fatalError(Here, '3 rotation angles must be given. Has only: '//numToChar(size(temp)))
      end if
      call self % setTransform(rotation=temp)
    end if

    ! Load radii and fill data
    call dict % get(radii, 'radii')
    call dict % get(fillNames, 'fills')

    ! Check values
    if (size(radii) /= size(fillNames)) then
      call fatalError(Here, 'Size of radii and fills does not match')

    else if (any(radii < ZERO)) then
      call fatalError(Here, 'Found -ve value of radius.')

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
    allocate(self % theta(self % nAz))
    allocate(self % planes(self % nAz/2))
    allocate(self % normals(self % nAz/2,2))
    dTheta = TWO_PI / self % nAz
    theta0 = HALF * dTheta
    do i = 1, self % nAz
      self % theta(i) = theta0 
      theta0 = theta0 + dTheta
    end do

    ! Build the planes, rotated at angle theta from the line x = 0
    do i = 1, self % nAz/2
      self % normals(i,1) = -sin(self % theta(i)) 
      self % normals(i,2) = cos(self % theta(i))
      call self % planes(i) % build(id=1, &
              norm = [self % normals(i,1), self % normals(i,2), ZERO], offset=ZERO)
    end do

    ! Create fill array
    allocate(fill(self % nAz * N))
    !allocate(fill(N))

    !do i = 1, N
    !    fill(i) = charToFill(fillNames(i), mats, Here)
    !end do
    do i = 1, N
      do j = 1, self % nAz
        fill((i-1)*self % nAz + j) = charToFill(fillNames(i), mats, Here)
      end do
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
    integer(shortInt)                       :: aIdx, rIdx, pIdx

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

    ! Find azimuthal segment
    do aIdx = 1, self % nAz
      if (aIdx > self % nAz/2) then
        pIdx = aIdx - self % nAz/2
        planeDir = -ONE
      else
        pIdx = aIdx
        planeDir = ONE
      end if
      ! Surface tolerance multiplier determined by relative direction of particle
      ! and theta
      if (planeDir*self % normals(pIdx,1)*u(1) + planeDir*self % normals(pIdx,2)*u(2) >= ZERO) then
        mul = -ONE
      else
        mul = ONE
      end if
      if (theta < self % theta(aIdx) + mul * self % planes(pIdx) % surfTol() ) exit
    end do
    ! If exceeded the search, theta is <2pi but greater than the largest theta.
    ! Therefore, it lies in the negative theta portion of the first segment.
    if (aIdx == self % nAz + 1) aIdx = 1
    localID = aIdx + self % nAz * (rIdx - 1)

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
    integer(shortInt)                     :: id, aIdx, rIdx, searchIdxP, searchIdxM
    integer(shortInt)                     :: sense, minus_sense, plus_sense
    character(100), parameter :: Here = 'distance (azimPinUniverse_class.f90)'

    ! Get local id
    id = coords % localID

    if (id < 1 .or. id > (size(self % r_sq) + 1) * self % nAz) then
      call fatalError(Here, 'Invalid local ID: '//numToChar(id))
    end if

    ! Identify annulus index and azimuthal index
    rIdx = floor((real(id-1)+0.1) / self % nAz) + 1
    aIdx = id - (self % nAz * (rIdx - 1))

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
    searchIdxP = aIdx
    if (aIdx > self % nAz/2) searchIdxP = searchIdxP - self % nAz/2

    ! Check for first or last cells circling round
    if (aIdx == 1) then
      minus_sense = MOVING_CLOCK_BACK
      plus_sense = MOVING_ANTI
    else if (aIdx == self % nAz) then
      plus_sense = MOVING_CLOCK_FORWARD
      minus_sense = MOVING_CLOCK
    else
      plus_sense = MOVING_ANTI
      minus_sense = MOVING_CLOCK
    end if

    ! If in the first cell or nAz/2 + 1 cell, find correct plane
    if (searchIdxP == 1) then
      searchIdxM = self % nAz/2
    else
      searchIdxM = searchIdxP - 1
    end if

    ! Identify which two planes (or only one if nAz = 2)
    ! Check to see if in second half of azimuthal segments
    if (self % nAz > 2) then
      d_plus  = self % planes(searchIdxP) % distance(coords % r, coords % dir)
      d_minus = self % planes(searchIdxM) % distance(coords % r, coords % dir)
    else
      d_plus  = self % planes(1) % distance(coords % r, coords % dir)
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
    character(100), parameter :: Here = 'cross (azimPinUniverse_class.f90)'

    ! Need to determine whether completing a circle, i.e., moving from
    ! the last azimuthal segment to the first and vice versa
    if (surfIdx == MOVING_CLOCK) then
      coords % localID = coords % localID - 1

    else if (surfIdx == MOVING_ANTI) then
      coords % localID = coords % localID + 1

    else if (surfIdx == MOVING_CLOCK_BACK) then
      coords % localID = coords % localID + self % nAz - 1

    else if (surfIdx == MOVING_CLOCK_FORWARD) then
      coords % localID = coords % localID - self % nAz + 1

    else if (surfIdx == MOVING_IN) then
      coords % localID = coords % localID - self % nAz

    else if (surfIdx == MOVING_OUT) then
      coords % localID = coords % localID + self % nAz

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
  elemental subroutine kill(self)
    class(azimPinUniverse), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Kill local
    if(allocated(self % r_sq)) deallocate(self % r_sq)
    if(allocated(self % annuli)) deallocate(self % annuli)
    if(allocated(self % theta)) deallocate(self % theta)
    if(allocated(self % planes)) deallocate(self % planes)
    if(allocated(self % normals)) deallocate(self % normals)
    self % nAz = 0

  end subroutine kill

end module azimPinUniverse_class
