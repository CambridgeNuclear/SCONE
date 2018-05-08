module coord_class
  use numPrecision
  use universalVariables
  use genericProcedures

  implicit none
  private

  type, public :: coord
    real(defReal), dimension(3) :: r                   ! position
    real(defReal), dimension(3) :: dir                 ! direction
    logical(defBool)            :: isRotated = .FALSE. ! is the co-ordinate in a rotated reference frame?
    integer(shortInt)           :: uniIdx = 0          ! index of the universe occupied
    integer(shortInt)           :: latIdx = 0          ! index of the lattice occupied
    integer(shortInt)           :: ijkIdx = 0          ! index of the lattice cell occupied (reduced to single digit)
    integer(shortInt)           :: cellIdx = 0         ! point to the cell occupied
  end type coord

  type, public :: coordList
    integer(shortInt)                          :: nesting ! depth of co-ordinate nesting
    type(coord), dimension(hardcoded_max_nest) :: lvl     ! array of coords nested successively deeper
  contains
    procedure :: init
    procedure :: addLevel
    procedure :: resetNesting
    procedure :: moveGlobal
    procedure :: moveLocal
    procedure :: rotate
    procedure :: cell
    procedure :: assignPosition
    procedure :: assignDirection
  end type coordList

contains

  !!
  !! Initialise co-ordinates with a global position
  !!
  subroutine init(self,r,dir)
    class(coordList), intent(inout)         :: self
    real(defReal), dimension(3), intent(in) :: r, dir
    self % lvl(1) % r = r
    self % lvl(1) % dir = dir
    self % nesting = 1
  end subroutine init

  !!
  !! Add another level of co-ordinates
  !!
  subroutine addLevel(self, offset, uniIdx, latIdx, ijkIdx)
    class(coordList), intent(inout)         :: self
    real(defReal), dimension(3), intent(in) :: offset
    integer(shortInt), intent(in)           :: uniIdx
    integer(shortInt), intent(in), optional :: latIdx
    integer(shortInt), intent(in), optional :: ijkIdx
    integer(shortInt)                       :: n

    n = self % nesting + 1
    self % lvl(n) % r = self % lvl(n-1) % r - offset
    self % lvl(n) % dir = self % lvl(n-1) % dir
    self % lvl(n) % uniIdx = uniIdx
    self % nesting = n
    if(present(latIdx)) then
      self % lvl(n) % latIdx = latIdx
      self % lvl(n) % ijkIdx = ijkIdx
    end if

  end subroutine addLevel

  !!
  !! Return the co-ordinates to the chosen level - if not specified, return to global
  !!
  subroutine resetNesting(self,n)
    class(coordList), intent(inout)         :: self
    integer(shortInt), intent(in), optional :: n
    integer(shortInt)                       :: nMin,nMax, i
    if(present(n)) then
      nMin = n
    else
      nMin = 1
    end if
    nMax = self % nesting
    do i = nMin+1,nMax
      self % lvl(i) % uniIdx  = 0
      self % lvl(i) % latIdx  = 0
      self % lvl(i) % cellIdx = 0
      self % lvl(i) % ijkIdx = 0
      self % lvl(i) % isRotated = .FALSE.
    end do
    self % nesting = nMin
  end subroutine resetNesting

  !!
  !! Move a point in global co-ordinates
  !!
  subroutine moveGlobal(self, distance)
    class(coordList), intent(inout) :: self
    real(defReal), intent(in)       :: distance
    call self % resetNesting()
    self % lvl(1) % r = self % lvl(1) % r + distance * self % lvl(1) % dir
  end subroutine moveGlobal

  !!
  !! Move a point in local co-ordinates down to nesting level n
  !!
  subroutine moveLocal(self, distance, n)
    class(coordList), intent(inout) :: self
    real(defReal), intent(in)       :: distance
    integer(shortInt), intent(in)   :: n
    integer(shortInt)               :: i
    call self % resetNesting(n)
    do i=1,n
      self % lvl(i) % r = self % lvl(i) % r + distance * self % lvl(i) % dir
    end do
  end subroutine moveLocal

  !!
  !! Rotate neutron direction
  !! Applies rotation vector to lower levels only if they are in a rotated geometry
  !! Otherwise, copies direction from the level above
  !!
  subroutine rotate(self,mu,phi)
    class(coordList), intent(inout) :: self
    real(defReal), intent(in)       :: mu
    real(defReal), intent(in)       :: phi
    integer(shortInt)               :: i

    ! Rotate directions in all nesting levels
    self % lvl(1) % dir = rotateVector(self % lvl(1) % dir, mu, phi)
    if (self % nesting > 1) then
      do i = 2,self % nesting
        if (self % lvl(i) % isRotated) then
          self % lvl(i) % dir = rotateVector(self % lvl(i) % dir, mu, phi)
        else
          self % lvl(i) % dir = self % lvl(i-1) % dir
        end if
      end do
    end if

  end subroutine rotate

  !!
  !! Returns the index of the cell occupied at the lowest level
  !!
  function cell(self)result(cellIdx)
    class(coordList), intent(in) :: self
    integer(shortInt)            :: cellIdx
    cellIdx = self % lvl(self % nesting) % cellIdx
  end function cell

  !!
  !! Assign the global position to an arbitrary value
  !!
  subroutine assignPosition(self, r)
    class(coordList), intent(inout)         :: self
    real(defReal), dimension(3), intent(in) :: r
    call self % resetNesting()
    self % lvl(1) % r = r
  end subroutine assignPosition

  !!
  !! Assign the global direction to an arbitrary value
  !!
  subroutine assignDirection(self, dir)
    class(coordList), intent(inout)         :: self
    real(defReal), dimension(3), intent(in) :: dir
    call self % resetNesting()
    self % lvl(1) % dir = dir
  end subroutine assignDirection


end module coord_class
