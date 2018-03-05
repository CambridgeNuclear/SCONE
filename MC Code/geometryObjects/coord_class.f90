module coord_class

  use numPrecision
  !use universalVariables
  use genericProcedures

  implicit none
  private

  !*** My version PARAMETER
  integer(shortInt), parameter :: max_nest = 3

  type, public :: coord
    real(defReal), dimension(3) :: r              ! position
    real(defReal), dimension(3) :: dir            ! direction
    integer(shortInt) :: uniInd = 0               ! index of the universe occupied
    integer(shortInt) :: latInd = 0               ! index of the lattice occupied
    !integer(shortInt) :: cellInd = 0              ! point to the cell occupied
  end type coord

  type, public :: coordList
    integer(shortInt) :: nesting               ! depth of co-ordinate nesting
    type(coord), dimension(max_nest) :: lvl    ! array of coords nested successively deeper
  contains
    procedure :: init
    procedure :: addLevel
    procedure :: resetNesting
    procedure :: moveGlobal
    procedure :: moveLocal
    procedure :: rotate
    procedure :: assignPosition
    procedure :: assignDirection
  end type coordList

contains

  !
  ! Initialise co-ordinates with a global position
  !
  subroutine init(self,r,dir)
    class(coordList), intent(inout) :: self
    real(defReal), dimension(3), intent(in) :: r, dir
    self % lvl(1) % r = r
    self % lvl(1) % dir = dir
    self % nesting = 1
  end subroutine init

  !
  ! Add another level of co-ordinates
  !
  subroutine addLevel(self, offset, uniInd, latInd)
    class(coordList), intent(inout)         :: self
    real(defReal), dimension(3), intent(in) :: offset
    integer(shortInt), intent(in)           :: uniInd
    integer(shortInt), intent(in),optional  :: latInd
    integer(shortInt)                       :: n
    n = self % nesting + 1
    self % lvl(n) % r = self % lvl(n-1) % r - offset
    self % lvl(n) % dir = self % lvl(n-1) % dir
    self % lvl(n) % uniInd = uniInd
    self % nesting = n
    if(present(latInd)) self % lvl(n) % latInd = latInd
  end subroutine addLevel

  !
  ! Add another level of co-ordinates on entering a lattice universe
  !
  !subroutine addLevelLattice(self, lat)
  !  class(coordList), intent(inout) :: self
  !  class(lattice_ptr), intent(in) :: lat
  !  real(defReal), dimension(3) :: r, u
  !  integer(shortInt) :: n
  !  integer(shortInt), dimension(3) :: ijkLat
  !  n = self % nesting + 1
  !  self % lvl(n) % lattice = lat
  !  r = self % lvl(n-1) % r
  !  u = self % lvl(n-1) % dir
  !  ijkLat = lat % findUniverse(r,u)
  !  self % lvl(n) % universe = lat % universes(ijkLat)
  !  self % lvl(n) % r = r - lat % localCoords(ijkLat) - self % lvl(n) % universe % offset()
  !  self % lvl(n) % dir = u
  !  self % nesting = n
  !end subroutine addLevelLattice

  !
  ! Return the co-ordinates to the chosen level - if not specified, return to global
  !
  subroutine resetNesting(self,n)
    class(coordList), intent(inout) :: self
    integer(shortInt), intent(in), optional :: n
    integer(shortInt) :: nMin,nMax, i

    if(present(n)) then
      nMin = n
    else
      nMin = 1
    end if
    nMax = self % nesting
    do i = nMin+1,nMax
      self % lvl(i) % uniInd = 0
      self % lvl(i) % latInd = 0
    end do
    self % nesting = nMin
  end subroutine resetNesting

  !!
  !! Move a point in global co-ordinates
  !!
  subroutine moveGlobal(self, distance)
    class(coordList), intent(inout) :: self
    real(defReal), intent(in) :: distance

    call self % resetNesting()
    self % lvl(1) % r = self % lvl(1) % r + distance * self % lvl(1) % dir
  end subroutine moveGlobal

  !!
  !! Move a point in local co-ordinates down to nesting level n
  !!
  subroutine moveLocal(self, distance, n)
    class(coordList), intent(inout) :: self
    real(defReal), intent(in) :: distance
    integer(shortInt), intent(in) :: n
    integer(shortInt) :: i

    call self % resetNesting(n)
    do i=1,n
      self % lvl(i) % r = self % lvl(i) % r + distance * self % lvl(i) % dir
    end do
  end subroutine moveLocal

  !!
  !! Rotate neutron direction
  !! Inefficient implementation, which does not account for not-rotated nesting levels
  !!
  subroutine rotate(self,mu,phi)
    class(coordList), intent(inout) :: self
    real(defReal), intent(in)       :: mu
    real(defReal), intent(in)       :: phi
    integer(shortInt)               :: i

    ! Rotate directions in all nesting levels
    do i =1,self % nesting
      self % lvl(i) % dir = rotateVector(self % lvl(i) % dir, mu, phi)

    end do

  end subroutine rotate


  !!
  !! Assign the global position to an arbitrary value
  !!
  subroutine assignPosition(self, r)
    class(coordList), intent(inout) :: self
    real(defReal), dimension(3), intent(in) :: r
    call self % resetNesting()
    self % lvl(1) % r = r
  end subroutine assignPosition

  !!
  !! Assign the global direction to an arbitrary value
  !!
  subroutine assignDirection(self, dir)
    class(coordList), intent(inout) :: self
    real(defReal), dimension(3), intent(in) :: dir

    call self % resetNesting()
    self % lvl(1) % dir = dir

  end subroutine assignDirection

end module coord_class
