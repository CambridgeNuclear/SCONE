module coord_class

  use numPrecision
  use universalVariables, only : HARDCODED_MAX_NEST, FP_REL_TOL
  use genericProcedures,  only : rotateVector, fatalError, numToChar

  implicit none
  private

  !!
  !! Co-ordinates in a single geometry level
  !!
  !! Co-ordinates are considered valid if:
  !!   * dir is normalised to 1.0 (norm2(dir) ~= 1.0)
  !!   * uniIdx, uniRootId & localID are set to +ve values
  !!
  !! Public Members:
  !!   r -> Position
  !!   dir -> Direction
  !!   isRotated -> Is rotated wrt previous (higher by 1) level
  !!   rotMat    -> Rotation matrix wrt previous level
  !!   uniIdx    -> Index of the occupied universe
  !!   uniRootID -> Location of the occupied universe in geometry graph
  !!   localID   -> Local cell in the occupied universe
  !!   cellIdx   -> Index of the occupied cell in cellShelf. 0 is cell is local to the universe
  !!
  !! Interface:
  !!   isValid -> True if co-ordinates are valid
  !!   display -> Print co-ordinates to the console
  !!   kill    -> Return to uninitialised state
  !!
  type, public :: coord
    real(defReal), dimension(3)   :: r         = ZERO
    real(defReal), dimension(3)   :: dir       = ZERO
    logical(defBool)              :: isRotated = .false.
    real(defReal), dimension(3,3) :: rotMat    = ZERO
    integer(shortInt)             :: uniIdx    = 0
    integer(shortInt)             :: uniRootID = 0
    integer(shortInt)             :: localID   = 0
    integer(shortInt)             :: cellIdx   = 0
  contains
    procedure :: isValid => isValid_coord
    procedure :: display => display_coord
    procedure :: kill    => kill_coord
  end type coord

  !!
  !! List of co-ordinates at diffrent level of a geometry
  !!
  !! Specifies the position of a particle in space
  !!
  !! It can exist in the following states:
  !!  ABOVE GEOMETRY     -> Nesting = 1. matIdx & uniqueID are < 0 (unassigned). Co-ordinates at
  !!                          level 1 are reliable.
  !!  PLACED IN GEOMETRY -> Nesting >=1. matIdx is assigned. Coordinates up to nesting are realible
  !!  UNINITIALISED      -> Is neither PLACED nor ABOVE
  !!
  !!  NOTE:
  !!   moveGlobal resets regionID & matIdx to 0
  !!   moveLocal  leaves regionID & matIdx unchanged
  !!
  !! Public Members:
  !!   nesting  -> Maximum currently occupied nexting level
  !!   lvl      -> Coordinates at each level
  !!   matIdx   -> Material index at the current position
  !!   uniqueID -> Unique cell ID at the current position
  !!
  !! Interface:
  !!   init              -> Initialise and place ABOVE GEOMETRY, given position and
  !!                          normalised direction
  !!   kill              -> Return to uninitialised state
  !!   isPlaced          -> True if co-ordinates are PLACED IN GEOMETRY
  !!   isAbove           -> True of co-ordinates are ABOVE GEOMETRY
  !!   isUninitialised   -> True if co-ordinates are UNINITIALISED
  !!   addLevel          -> Increment number of occupied levels by 1.
  !!   decreaseLevel     -> Decrease nesting to a lower level n.
  !!   takeAboveGeometry -> Change state to ABOVE GEOMETRY
  !!   moveGlobal        -> Move point along direction ABOVE GEOMETRY
  !!   moveLocal         -> Move point along direction above and including level n
  !!   rotate            -> Rotate by the cosine of polar deflection mu, and azimuthal angle phi
  !!   cell              -> Return cellIdx at the lowest level
  !!   assignPosition    -> Set position and take ABOVE GEOMETRY
  !!   assignDirection   -> Set direction and do not change coordList state
  !!
  type, public :: coordList
    integer(shortInt)                          :: nesting = 0
    type(coord), dimension(HARDCODED_MAX_NEST) :: lvl
    integer(shortInt)                          :: matIdx   = -3
    integer(shortInt)                          :: uniqueId = -3
  contains
    ! Build procedures
    procedure :: init
    procedure :: kill

    ! State enquiry procedures
    procedure :: isPlaced
    procedure :: isAbove
    procedure :: isUninitialised

    ! Interface procedures
    procedure :: addLevel
    procedure :: decreaseLevel
    procedure :: takeAboveGeom
    procedure :: moveGlobal
    procedure :: moveLocal
    procedure :: rotate
    procedure :: cell
    procedure :: assignPosition
    procedure :: assignDirection

  end type coordList

contains

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! coord Procedures
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Returns .true. if coordinates are valid
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   True if coord is valid. See type doc-comment for definition of valid.
  !!
  elemental function isValid_coord(self) result(correct)
    class(coord), intent(in) :: self
    logical(defBool)         :: correct

    ! Direction vector is normalised within floating point tolerance
    correct = abs(norm2(self % dir) - ONE) < FP_REL_TOL

    correct = correct .and. self % uniIdx  > 0
    correct = correct .and. self % localID > 0
    correct = correct .and. self % uniRootId > 0

  end function isValid_coord

  !!
  !! Print to screen contents of the coord
  !!
  subroutine display_coord(self)
    class(coord), intent(in) :: self

    print *, "R: ", self % r
    print *, "U: ", self % dir
    print *, "UniIdx: ", numToChar(self % uniIDx), " LocalID: ", numToChar(self % localID), &
             "UniRootId", numToChar(self % uniRootID)

  end subroutine display_coord

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill_coord(self)
    class(coord), intent(inout) :: self

    self % r         = ZERO
    self % dir       = ZERO
    self % isRotated = .false.
    self % rotMat    = ZERO
    self % uniIdx    = 0
    self % uniRootID = 0
    self % localID   = 0
    self % cellIdx   = 0

  end subroutine kill_coord


!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! coordList Procedures
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Initialise coordList
  !!
  !! Change state from UNINITIALISED to ABOVE GEOMETRY
  !!
  !! Args:
  !!   r [in] -> Position in level 1
  !!   u [in] -> Normalised direction in level 1 (norm2(u)=1.0)
  !!
  !! NOTE:
  !!   Does not check if u is normalised!
  !!
  pure subroutine init(self, r, u)
    class(coordList), intent(inout)         :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u

    call self % takeAboveGeom()

    self % lvl(1) % r = r
    self % lvl(1) % dir = u
    self % nesting = 1

  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(coordList), intent(inout) :: self

    self % nesting  = 0
    self % matIdx   = -3
    self % uniqueID = -3

    ! Kill coordinates
    call self % lvl % kill()

  end subroutine kill

  !!
  !! Return true if co-ordinates List is placed in geometry
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   True if co-ordinates are PLACED.
  !!
  elemental function isPlaced(self) result(isIt)
    class(coordList), intent(in) :: self
    logical(defBool)             :: isIt

    isIt = (self % matIdx > 0) .and. (self % uniqueID > 0) .and. (self % nesting >= 1)

  end function isPlaced

  !!
  !! Return true if co-ordinates are above geometry
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   True if co-ordinates are ABOVE GEOMETRY
  !!
  elemental function isAbove(self) result(isIt)
    class(coordList), intent(in) :: self
    logical(defBool)             :: isIt

    isIt = (self % matIdx < 0) .and. (self % uniqueID < 0) .and. (self % nesting == 1)

  end function isAbove

  !!
  !! Return true if coordinates are uninitialised
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   True if co-ordinates are UNINITIALISED
  !!
  elemental function isUninitialised(self)  result(isIt)
    class(coordList), intent(in) :: self
    logical(defBool)             :: isIt

    isIt = .not.( self % isPlaced() .or. self % isAbove() )

  end function isUninitialised

  !!
  !! Add another level of co-ordinates
  !!
  !! Simply increments nesting counter
  !!
  !! Args:
  !!   None
  !!
  pure subroutine addLevel(self)
    class(coordList), intent(inout) :: self

    self % nesting = self % nesting + 1

  end subroutine addLevel

  !!
  !! Takes co-ordinates above the geometry
  !!
  !! State changes to ABOVE GEOMETRY
  !!
  !! Args:
  !!   None
  !!
  !! NOTE:
  !!   If called on UNINITIALISED may result in unnormalised direction at level 1!
  !!
  elemental subroutine takeAboveGeom(self)
    class(coordList), intent(inout) :: self

    self % nesting  = 1
    self % matIdx   = -3
    self % uniqueID = -3

  end subroutine takeAboveGeom

  !!
  !! Decrease nestting to level n
  !!
  !! Args:
  !!   n [in] -> New nesting level
  !!
  !! Errors:
  !!   fatalError if n is -ve or larger then current nesting
  !!
  subroutine decreaseLevel(self, n)
    class(coordList), intent(inout)  :: self
    integer(shortInt), intent(in)    :: n
    character(100),parameter :: Here = 'decreaseLevel (coord_class.f90)'

    if (n > self % nesting .or. n < 1) then
      call fatalError(Here,'New nesting: '//numToChar(n)//' is invalid. Current nesting: '//&
                            numToChar(self % nesting))
    end if

    self % nesting = n

  end subroutine decreaseLevel

  !!
  !! Move a point ABOVE the geometry
  !!
  !! Changes state to ABOVE GEOMETRY
  !!
  !! Args:
  !!   d [in] -> Distance (+ve or -ve)
  !!
  !! Errors:
  !!   If d < 0 then movment is backwards.
  !!
  elemental subroutine moveGlobal(self, d)
    class(coordList), intent(inout) :: self
    real(defReal), intent(in)       :: d

    call self % takeAboveGeom()
    self % lvl(1) % r = self % lvl(1) % r + d * self % lvl(1) % dir

  end subroutine moveGlobal

  !!
  !! Move point inside the geometry
  !!
  !! Moves above and including level n
  !! Does not change matIdx nor uniqueID
  !!
  !! Args:
  !!   d [in] -> Distance (+ve or -ve)
  !!   n [in] -> Nesting level
  !!
  !! Errors:
  !!   If d < 0.0 movment is backwards
  !!
  subroutine moveLocal(self, d, n)
    class(coordList), intent(inout) :: self
    real(defReal), intent(in)       :: d
    integer(shortInt), intent(in)   :: n
    integer(shortInt)               :: i

    call self % decreaseLevel(n)
    do i = 1 , n
      self % lvl(i) % r = self % lvl(i) % r + d * self % lvl(i) % dir
    end do

  end subroutine moveLocal

  !!
  !! Rotate direction of the point
  !!
  !! Does not change the state of co-ordinates
  !!
  !! Args:
  !!   mu [in]  -> Cosine of polar deflection angle <-1,1>
  !!   phi [in] -> Azimuthal deflection angle <0;2*pi>
  !!
  elemental subroutine rotate(self, mu, phi)
    class(coordList), intent(inout) :: self
    real(defReal), intent(in)       :: mu
    real(defReal), intent(in)       :: phi
    integer(shortInt)               :: i

    ! Rotate directions in all nesting levels
    self % lvl(1) % dir = rotateVector(self % lvl(1) % dir, mu, phi)

    ! Propagate rotation to lower levels
    do i = 2, self % nesting
      if (self % lvl(i) % isRotated) then
        ! Note that rotation must be performed with the matrix
        ! Deflections by mu & phi depend on coordinates
        ! Deflection by the same my & phi may be diffrent at diffrent, rotated levels!
        self % lvl(i) % dir = matmul(self % lvl(i) % rotMat, self % lvl(i-1) % dir)

      else
        self % lvl(i) % dir = self % lvl(i-1) % dir

      end if
    end do

  end subroutine rotate

  !!
  !! Returns the index of the cell occupied at the lowest level
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   cellIdx at the lowest ocupied level
  !!
  elemental function cell(self)result(cellIdx)
    class(coordList), intent(in) :: self
    integer(shortInt)            :: cellIdx

    cellIdx = self % lvl(max(self % nesting, 1)) % cellIdx

  end function cell

  !!
  !! Take co-ordinates ABOVE GEOMETRY and assign new position
  !!
  !! Args:
  !!   r [in] -> New position at level 1
  !!
  pure subroutine assignPosition(self, r)
    class(coordList), intent(inout)         :: self
    real(defReal), dimension(3), intent(in) :: r

    call self % takeAboveGeom()
    self % lvl(1) % r = r

  end subroutine assignPosition

  !!
  !! Assign new direction
  !!
  !! Does not change the state of co-ordinates
  !!
  !! Args:
  !!   u [in] -> New normalised direction at level 1
  !!
  !! NOTE:
  !!   Does not check if u is normalised!
  !!
  pure subroutine assignDirection(self, u)
    class(coordList), intent(inout)         :: self
    real(defReal), dimension(3), intent(in) :: u
    integer(shortInt)                       :: i

    ! Assign new direction in global frame
    self % lvl(1) % dir = u

    ! Propage the change to lower levels
    do i = 2, self % nesting
      if(self % lvl(i) % isRotated) then
        self % lvl(i) % dir = matmul(self % lvl(i) % rotMat, self % lvl(i-1) % dir)

      else
        self % lvl(i) % dir = self % lvl(i-1) % dir

      end if
    end do

  end subroutine assignDirection

end module coord_class
