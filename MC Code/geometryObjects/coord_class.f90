module coord_class

  use numPrecision
!  use universalVariables
  use genericProcedures

 ! use cell_class
 ! use universe_class
 ! use lattice_class

  implicit none
  private

  integer(shortInt), parameter :: max_nest = 5

  type, public :: coord
    real(defReal), dimension(3) :: r              ! position
    real(defReal), dimension(3) :: dir            ! direction
  !  type(cell_ptr) :: cell                        ! point to the cell occupied
  !  type(universe_ptr) :: universe                ! point to the universe occupied
  !  type(lattice_ptr) :: lattice                  ! point to the lattice occupied
  !  logical(defBool) :: inLattice                 ! whether point is inside a lattice cell
  end type coord

  type, public :: coordList
    integer(shortInt) :: nesting               ! depth of co-ordinate nesting
    type(coord), dimension(max_nest) :: lvl    ! array of coords nested successively deeper
  contains
    procedure :: init
   ! procedure :: addLevelUniverse
    !procedure :: addLevelLattice
    !procedure :: resetNesting
  end type coordList

contains

  !
  ! Initialise co-ordinates with a global position
  !
  subroutine init(self,r,dir)
    implicit none
    class(coordList), intent(inout) :: self
    real(defReal), dimension(3), intent(in) :: r, dir
    self % lvl(1) % r = r
    self % lvl(1) % dir = dir
    self % nesting = 1
  end subroutine init

  !
  ! Add another level of co-ordinates offset on entering a universe
  !
!  subroutine addLevelUniverse(self, uni)
!    implicit none
!    class(coordList), intent(inout) :: self
!    class(universe_ptr), intent(in) :: uni
!    integer(shortInt) :: n
!    n = self % nesting + 1
!    self % lvl(n) % r = self % lvl(n-1) % r - uni % offset()
!    self % lvl(n) % dir = self % lvl(n-1) % dir
!    self % lvl(n) % universe = uni
!    self % nesting = n
!  end subroutine addLevelUniverse

  !
  ! Add another level of co-ordinates on entering a lattice universe
  !
!  subroutine addLevelLattice(self, lat)
!    implicit none
!    class(coordList), intent(inout) :: self
!    class(lattice_ptr), intent(in) :: lat
!    real(defReal), dimension(3) :: r, u
!    integer(shortInt) :: n
!    integer(shortInt), dimension(3) :: ijkLat
!    n = self % nesting + 1
!    self % lvl(n) % lattice = lat
!    r = self % lvl(n-1) % r - lat % offset()
!    u = self % lvl(n-1) % dir
!    ijkLat = lat % findUniverse(r,u)
!    self % lvl(n) % universe = lat % universes(ijkLat)
!    self % lvl(n) % r = r - lat % localCoords(ijkLat) - self % lvl(n) % universe % offset()
!    self % lvl(n) % dir = u
!    self % lvl(n) % inLattice = .true.
!    self % nesting = n
!  end subroutine addLevelLattice

  !
  ! Return the co-ordinates to the chosen level - if not specified, return to global
  !
!  subroutine resetNesting(self,n)
!    implicit none
!    class(coordList), intent(inout) :: self
!    integer(shortInt), intent(in), optional :: n
!    integer(shortInt) :: nMin,nMax, i
!
!    if(present(n)) then
!      nMin = n
!    else
!      nMin = 1
!    end if
!
!    nMax = self % nesting
!    do i = nMin+1,nMax
!      call self % lvl(i) % cell % nullify()
!      call self % lvl(i) % universe % nullify()
!      call self % lvl(i) % lattice % nullify()
!      self % lvl(i) % inLattice = .false.
!    end do
!    self % nesting = nMin
!  end subroutine resetNesting

end module coord_class
