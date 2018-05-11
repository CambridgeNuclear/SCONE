!
! The lattice class contains an array of universes in a cubic array.
!
! A particle can be quickly and efficiently located within a sub-universe through
! co-ordinate comparison and indexing.
!
module lattice_class
  use numPrecision
  use genericProcedures
  use universalVariables
  use cell_class
  use universe_class

  implicit none
  private

  type, public :: lattice
    type(universe_ptr), dimension(:,:,:), allocatable :: universes         ! universes contained by lattice in 3D
    real(defReal), dimension(3)                       :: pitch             ! pitch of universes in each dimension
    real(defReal), dimension(3)                       :: corner            ! bottom left corner of the lattice in the parent cell's universe
    integer(shortInt), dimension(3)                   :: extent            ! number of cells in each direction
    integer(shortInt)                                 :: totalCells        ! total number of cells in the lattice
    integer(shortInt)                                 :: id                ! unique ID identifying the lattice
    logical(defBool)                                  :: is3D              ! simplifies lattice procedures if it is not 3D
    integer(shortInt)                                 :: outsideMatIdx = 0 ! material contained in any index outside the lattice
    character(100), public                            :: name = ""
  contains
    procedure :: init
    procedure :: findUniverse  ! return the indices of the daughter universe given a point
    procedure :: localCoords   ! return the local co-ordinates of a given universe
    procedure :: getDistance   ! assuming a point is within a lattice cell, find the distance to the boundary
    procedure :: insideLattice ! checks that an index is within the lattice
    procedure :: getijkIdx     ! given a lattice array position, calculate a unique index given lattice dimensions
  end type lattice

  type, public :: lattice_ptr
    class(lattice), pointer :: ptr
  contains
    procedure :: init => init_ptr
    procedure :: findUniverse => findUniverse_ptr   ! return the indices of the daughter universe given a point
    procedure :: localCoords => localCoords_ptr     ! return the local co-ordinates of a given universe
    procedure :: getDistance => getDistance_ptr     ! assuming a point is within a lattice cell, find the distance to the boundary
    procedure :: insideLattice => insideLattice_ptr ! ensures that a given index is within the lattice
    procedure :: getijkIdx => getijkIdx_ptr         ! get unique ijkIdx given lattice array position
    procedure :: universes => universes_ptr         ! returns a universe of the lattice which is pointed to if given an index
    procedure :: name => name_ptr
    !procedure :: outsideMatInd => outsideMatInd_ptr ! returns the material index of a region outside the lattice
    procedure :: kill
    generic   :: assignment(=) => lattice_ptr_assignment, lattice_ptr_assignment_target
    procedure,private :: lattice_ptr_assignment
    procedure,private :: lattice_ptr_assignment_target
  end type lattice_ptr

contains

  !!
  !! Initialise the lattice by giving its dimensions and the
  !! array of universes which it contains
  !!
  !! Must also create the cell which will contain each of the universes
  !! These universes are then assigned this cell as their parent cell
  !!
  subroutine init(self, pitch, corner, universes, id, is3D, name)
    class(lattice), intent(inout)                     :: self
    real(defReal), dimension(3), intent(in)           :: pitch
    real(defReal), dimension(3), intent(in)           :: corner
    class(universe_ptr), dimension(:,:,:), intent(in) :: universes
    integer(shortInt), intent(in)                     :: id
    character(*), intent(in), optional                :: name
    logical(defBool), intent(in)                      :: is3D
    integer(shortInt), dimension(3)                   :: sz

    sz(1) = size(universes,1)
    sz(2) = size(universes,2)
    sz(3) = size(universes,3)
    self % pitch = pitch
    self % extent = sz
    self % corner = corner
    self % is3D = is3D
    self % id = id
    if(is3D) then
      if(any(pitch < surface_tol)) &
      call fatalError('init, lattice','3D lattice pitches must be greater than surface tolerance')
      allocate(self % universes(sz(1),sz(2),sz(3)))
      self % totalCells = sz(1) * sz(2) * sz(3)
    else
      if((pitch(1) < surface_tol) .OR. (pitch(2) < surface_tol)) &
      call fatalError('init, lattice','x and y lattice pitches must be greater than surface tolerance')
      allocate(self % universes(sz(1),sz(2),1))
      self % universes = universes(:,:,1:1)
      self % totalCells = sz(1) * sz(2)
    end if

    if(present(name)) self % name = name
  end subroutine init

  !!
  !! Assuming coords are in the lattice system, find the indices of the unvierse the point occupies
  !!
  function findUniverse(self,r,u) result(ijk)
    class(lattice), intent(in)              :: self
    real(defReal), intent(in), dimension(3) :: r
    real(defReal), intent(in), dimension(3) :: u
    real(defReal), dimension(3)             :: corner, pitch
    integer(shortInt), dimension(3)         :: ijk

    corner = self % corner
    pitch = self % pitch

    ijk(1) = ceiling((r(1) - corner(1))/pitch(1))
    ijk(2) = ceiling((r(2) - corner(2))/pitch(2))
    if (self % is3D) then
      ijk(3) = ceiling((r(3) - corner(3))/pitch(3))
    else
      ijk(3) = 1
    end if

  end function findUniverse

  !!
  !! Returns the co-ordinate offset for a particular universe in the lattice
  !!
  function localCoords(self, ijk) result(shift)
    class(lattice), intent(in)                  :: self
    integer(shortInt), dimension(3), intent(in) :: ijk
    real(defReal), dimension(3)                 :: shift

    shift(1) = self % corner(1) + (ijk(1) - HALF)*self % pitch(1)
    shift(2) = self % corner(2) + (ijk(2) - HALF)*self % pitch(2)
    if (self % is3D) then
      shift(3) = self % corner(3) + (ijk(3) - HALF)*self % pitch(3)
    else
      shift(3) = ZERO
    end if

  end function localCoords

  !!
  !! Returns the distance that a particle must travel to the boundary of a lattice cell
  !! Assumes the particle is in local co-ordinates
  !!
  function getDistance(self, r, u) result(distance)
    class(lattice), intent(in)              :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    real(defReal), dimension(3)             :: bound
    real(defReal)                           :: distance
    real(defReal)                           :: testDistance
    integer(shortInt)                       :: i, maxDim

    distance = INFINITY
    ! Find the bounds in each direction which the particle must intersect
    bound(1) = sign(self % pitch(1)*HALF, u(1))
    bound(2) = sign(self % pitch(2)*HALF, u(2))
    if (self % is3D) then
      bound(3) = sign(self % pitch(3)*HALF, u(3))
      maxDim = 3
    else
      maxDim = 2
    end if

    do i = 1,maxDim
      if(abs(r(i) - bound(i)) < surface_tol) then
        testDistance = INFINITY
      elseif (u(i) == ZERO) then
        testDistance = INFINITY
      else
        testDistance = (bound(i) - r(i))/u(i)
      end if
      if (distance > testDistance) distance = testDistance
    end do

  end function getDistance

  !!
  !! Ensures that the index provided is actually inside the lattice
  !!
  function insideLattice(self,ijk) result(isInside)
    class(lattice), intent(in)                  :: self
    integer(shortInt), dimension(3), intent(in) :: ijk
    logical(defBool)                            :: isInside
    isInside = all((ijk>0).AND.(ijk<=self%extent))
  end function insideLattice

  !!
  !! Return a single digit identifier for a lattice position given its 3D index
  !!
  function getijkIdx(self,ijk) result(ijkIdx)
    class(lattice), intent(in)                  :: self
    integer(shortInt), dimension(3), intent(in) :: ijk
    integer(shortInt)                           :: ijkIdx

    ijkIdx = ijk(1) + self % extent(1) * (ijk(2) - 1 + self % extent(2) * (ijk(3) - 1))

  end function getijkIdx

!!
!! Pointer wrapper procedures
!!

  !!
  !! Initialise the lattice by giving its dimensions and the
  !! array of universes which it contains
  !!
  !! Must also create the cell which will contain each of the universes
  !! These universes are then assigned this cell as their parent cell
  !!
  subroutine init_ptr(self, pitch, corner, universes, id, is3D, name)
    class(lattice_ptr), intent(inout)                 :: self
    real(defReal), dimension(3), intent(in)           :: pitch
    real(defReal), dimension(3), intent(in)           :: corner
    class(universe_ptr), dimension(:,:,:), intent(in) :: universes
    integer(shortInt), intent(in)                     :: id
    logical(defBool), intent(in)                      :: is3D
    character(*), intent(in)                          :: name
    call self % ptr % init(pitch, corner, universes, id, is3D, name)
  end subroutine init_ptr

  !!
  !! Assuming coords are in the lattice system, find the indices of the unvierse the point occupies
  !!
  function findUniverse_ptr(self,r,u)result(ijk)
    class(lattice_ptr), intent(in)          :: self
    real(defReal), intent(in), dimension(3) :: r
    real(defReal), intent(in), dimension(3) :: u
    integer(shortInt), dimension(3)         :: ijk
    ijk = self % ptr % findUniverse(r,u)
  end function findUniverse_ptr

  !!
  !! Returns the co-ordinate offset for a particular universe in the lattice
  !!
  function localCoords_ptr(self, ijk)result(shift)
    class(lattice_ptr), intent(in)              :: self
    integer(shortInt), dimension(3), intent(in) :: ijk
    real(defReal), dimension(3)                 :: shift
    shift = self % ptr % localCoords(ijk)
  end function localCoords_ptr

  !!
  !! Returns the distance that a particle must travel to the boundary of a lattice cell
  !! Assumes the particle is in local co-ordinates
  !!
  function getDistance_ptr(self, r, u)result(distance)
    class(lattice_ptr), intent(in)          :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    real(defReal)                           :: distance
    distance = self % ptr % getDistance(r,u)
  end function getDistance_ptr

  !!
  !! Ensures that the index of the universe provided is within the lattice
  !!
  function insideLattice_ptr(self,ijk)result(isInside)
    class(lattice_ptr), intent(in)              :: self
    integer(shortInt), dimension(3), intent(in) :: ijk
    logical(defBool)                            :: isInside
    isInside = self % ptr % insideLattice(ijk)
  end function insideLattice_ptr

  !!
  !! Returns the universe given a specific lattice index
  !!
  function universes_ptr(self,ijk)result(uni)
    class(lattice_ptr), intent(in)              :: self
    integer(shortInt), dimension(3), intent(in) :: ijk
    type(universe_ptr)                          :: uni
    uni = self % ptr % universes(ijk(1),ijk(2),ijk(3))
  end function universes_ptr

  !!
  !! Return a single digit identifier for a lattice position given its 3D index
  !!
  function getijkIdx_ptr(self,ijk) result(ijkIdx)
    class(lattice_ptr), intent(in)              :: self
    integer(shortInt), dimension(3), intent(in) :: ijk
    integer(shortInt)                           :: ijkIdx
    ijkIdx = self % ptr % getijkIdx(ijk)
  end function getijkIdx_ptr

  !!
  !! Returns the name of the lattice
  !!
  function name_ptr(self)result(name)
    class(lattice_ptr), intent(in) :: self
    character(100)                 :: name
    name = self % ptr % name
  end function name_ptr

  !!
  !! Returns the material index of a point outside the lattice
  !!
  !function outsideMaterialInd_ptr(self)result(outsideMaterialInd)
  !  class(lattice_ptr), intent(in) :: self
  !  integer(shortInt) :: outsideMaterialInd
  !  outsideMaterialInd = self % ptr % outsideMaterialInd
  !  if (outsideMaterialInd ==0) &
  !  call fatalError('outsideMaterialInd_ptr, lattice','Lattice does not have an outside material')
  !end function outsideMaterialInd_ptr

  subroutine lattice_ptr_assignment(LHS,RHS)
    class(lattice_ptr), intent(out)  :: LHS
    type(lattice_ptr), intent(in)    :: RHS

    if(associated(LHS % ptr)) deallocate(LHS % ptr)
    LHS % ptr => RHS % ptr
  end subroutine lattice_ptr_assignment

  subroutine lattice_ptr_assignment_target(LHS,RHS)
    class(lattice_ptr), intent(out)        :: LHS
    class(lattice), target, intent(in)     :: RHS

    if(associated(LHS % ptr)) deallocate(LHS % ptr)
    LHS % ptr => RHS
  end subroutine lattice_ptr_assignment_target

  subroutine kill(self)
    class(lattice_ptr), intent(inout) :: self
    self % ptr => null()
  end subroutine kill

end module lattice_class
