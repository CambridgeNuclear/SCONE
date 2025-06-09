module pinUniverse_class

  use numPrecision
  use universalVariables, only : INF, targetNotFound
  use genericProcedures,  only : fatalError, numToChar, swap
  use dictionary_class,   only : dictionary
  use coord_class,        only : coord
  use charMap_class,      only : charMap
  use surfaceShelf_class, only : surfaceShelf
  use cylinder_class,     only : cylinder
  use cell_inter,         only : cell
  use cellShelf_class,    only : cellShelf
  use universe_inter,     only : universe, kill_super => kill, charToFill
  implicit none
  private

  ! Parameters
  ! Are public for use in unit tests
  integer(shortInt), parameter, public :: MOVING_IN = -1, MOVING_OUT = -2

  !!
  !! Universe that represents a single pin
  !!
  !! Is composed from co-centring cylinders. Central cell has local ID 1 and the ID
  !! increases with subsequent rings.
  !!
  !! Sample Dictionary Input:
  !!   pinUni {
  !!     id 7;
  !!     type pinUniverse;
  !!     #origin (1.0 0.0 0.1);#
  !!     #rotation (30.0 0.0 0.0);#
  !!     radii (3.0 4.5 0.0 1.0 );
  !!     fills (u<3> void clad u<4>);
  !!   }
  !!
  !!  There must be 0.0 entry, which indicates outermost annulus (infinite radius).
  !!  `fills` and `radii` are given as pairs by position in the input arrays. Thus, fills
  !!  are sorted together with the `radii`. As a result, in the example, local cell 1 is
  !!  filled with u<4>, cell 2 with u<3> etc.
  !!
  !! Public Members:
  !!  r_sqr  -> Array of radius^2 for each annulus
  !!  annuli -> Array of cylinder surfaces that represent different annuli
  !!
  !! Interface:
  !!   universe interface
  !!
  type, public, extends(universe) :: pinUniverse
    private
    real(defReal), dimension(:), allocatable  :: r_sq
    type(cylinder), dimension(:), allocatable :: annuli
  contains
    ! Superclass procedures
    procedure :: init
    procedure :: kill
    procedure :: findCell
    procedure :: distance
    procedure :: cross
    procedure :: cellOffset
    procedure :: getNormal
  end type pinUniverse

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
    class(pinUniverse), intent(inout)                        :: self
    integer(shortInt), dimension(:), allocatable, intent(out) :: fill
    class(dictionary), intent(in)                             :: dict
    type(cellShelf), intent(inout)                            :: cells
    type(surfaceShelf), intent(inout)                         :: surfs
    type(charMap), intent(in)                                 :: mats
    integer(shortInt)                             :: idx, N, i
    real(defReal), dimension(:), allocatable      :: radii
    character(nameLen), dimension(:), allocatable :: fillNames
    character(100), parameter :: Here = 'init (pinUniverse_class.f90)'

    ! Setup the base class
    ! With: id, origin rotations...
    call self % setupBase(dict)

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
    if (radii(idx) /= ZERO) call fatalError(Here, 'Did not found outermost element with radius 0.0.')
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

    ! Create fill array
    allocate(fill(N))

    do i = 1, N
      fill(i) = charToFill(fillNames(i), mats, Here)
    end do


  end subroutine init

  !!
  !! Find local cell ID given a point
  !!
  !! See universe_inter for details.
  !!
  subroutine findCell(self, localID, cellIdx, r, u)
    class(pinUniverse), intent(inout)       :: self
    integer(shortInt), intent(out)          :: localID
    integer(shortInt), intent(out)          :: cellIdx
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    real(defReal)                           :: r_sq, mul

    r_sq = r(1)*r(1) + r(2)*r(2)
    cellIdx = 0

    ! Need to include surface tolerance. Determine multiplier by direction
    if ( r(1)*u(1) + r(2)*u(2) >= ZERO) then
      mul = -ONE
    else
      mul = ONE
    end if

    ! Find local cell
    do localID = 1, size(self % r_sq)
      if( r_sq < self % r_sq(localID) + mul * self % annuli(localID) % surfTol() ) return

    end do
    ! If reached here localID = size(self % r_sq) + 1

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
    class(pinUniverse), intent(inout)  :: self
    real(defReal), intent(out)         :: d
    integer(shortInt), intent(out)     :: surfIdx
    type(coord), intent(in)            :: coords
    real(defReal)                      :: d_out, d_in
    integer(shortInt)                  :: id
    character(100), parameter :: Here = 'distance (pinUniverse_class.f90)'

    ! Get local id
    id = coords % localID

    if (id < 1 .or. id > size(self % r_sq) + 1) then
      call fatalError(Here, 'Invalid local ID: '//numToChar(id))
    end if

    ! Outer distance
    if (id > size(self % r_sq)) then
      d_out = INF

    else
      d_out = self % annuli(id) % distance(coords % r, coords % dir)
    end if

    ! Inner distance
    if (id == 1) then
      d_in = INF
    else
      d_in = self % annuli(id-1) % distance(coords % r, coords % dir)
    end if

    ! Select distance and surface
    if ( d_in < d_out) then
      surfIdx = MOVING_IN
      d = d_in

    else
      surfIdx = MOVING_OUT
      d = d_out
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
    class(pinUniverse), intent(inout) :: self
    type(coord), intent(inout)        :: coords
    integer(shortInt), intent(in)     :: surfIdx
    character(100), parameter :: Here = 'cross (pinUniverse_class.f90)'

    if (surfIdx == MOVING_IN) then
      coords % localID = coords % localID - 1

    else if (surfIdx == MOVING_OUT) then
      coords % localID = coords % localID + 1

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
    class(pinUniverse), intent(in) :: self
    type(coord), intent(in)         :: coords
    real(defReal), dimension(3)     :: offset

    ! There is no cell offset
    offset = ZERO

  end function cellOffset
  
  !!
  !! Return normal for the given surfIdx at a given point
  !!
  !! See universe_inter for details.
  !!
  function getNormal(self, surfIdx, coords) result (normal)
    class(pinUniverse), intent(in) :: self
    integer(shortInt), intent(in)  :: surfIdx
    type(coord), intent(in)        :: coords
    real(defReal), dimension(3)    :: normal
    integer(shortInt)               :: cIdx
    character(100), parameter :: Here = 'getNormal (pinUniverse_class.f90)'

    ! Local ID and surfIdx should be sufficient to identify which annulus
    ! is being crossed and therefore which normal to produce
    cIdx = coords % localID

    if (surfIdx == MOVING_IN) then
      cIdx = cIdx
    elseif (surfIdx == MOVING_OUT) then
      cIdx = cIdx - 1
    else
      call fatalError(Here, 'Unrecognised surfIdx: '//numToChar(surfIdx))
    end if

    if (cIdx > size(self % annuli) .or. cIdx < 1) call fatalError(Here,&
            'Invalid cIdx requested: '// numToChar(cIdx))

    normal = self % annuli(cIdx) % normal(coords % r, coords % dir)   

  end function getNormal

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(pinUniverse), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Kill local
    if(allocated(self % r_sq)) deallocate(self % r_sq)
    if(allocated(self % annuli)) deallocate(self % annuli)

  end subroutine kill

end module pinUniverse_class
