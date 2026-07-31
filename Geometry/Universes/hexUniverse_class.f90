module hexUniverse_class

  use numPrecision
  use universalVariables,     only : INF, SURF_TOL
  use genericProcedures,      only : fatalError, numToChar, swap
  use dictionary_class,       only : dictionary
  use coord_class,            only : coord
  use charMap_class,          only : charMap
  use surfaceShelf_class,     only : surfaceShelf
  use box_class,              only : box
  use cell_inter,             only : cell
  use cellShelf_class,        only : cellShelf
  use universe_inter,         only : universe, kill_super => kill, charToFill
  use hexagonalLattice_class, only : hexagonalLattice

  implicit none
  private

  ! Parameters
  ! Note FACE<n>_POS/NEG and AX_MIN/AX_MAX are provided by hexagonalLattice_class.f90
  ! Outline must be different from those (which go from -1 to -8)
  ! Public so it can be accessed by unit tests
  integer(shortInt), public, parameter :: OUTLINE_SURF = -9

  ! Options for the offset map
  integer(shortInt), parameter :: local = 1, noOffset = 0

  !!
  !! Hexagonal lattice, stacked in layers along a Cartesian axis
  !!
  !! Universe consists of a hexagonal arrangement of cells (e.g. 5 x 5), which may be
  !! further stacked into layers along the axis perpendicular to the hexagonal plane
  !! (e.g. giving a 5 x 5 x 2 lattice). Centre of the lattice is placed at the origin.
  !! An additional cell is placed beyond the lattice called background (or out) cell.
  !!
  !! Local ID is 1 at the (i,j,k) = (1,1,1) cell. It increases first with i (1st in-plane
  !! oblique co-ordinate), then j (2nd in-plane oblique co-ordinate) and lastly k (the axial
  !! layer). Cells inside the lattice can only be filled with a universe (given as integer
  !! ID). Background cell can have any filling given by keyword (material or universe).
  !!
  !! Every lattice cell has an offset to its centre (so the centre of the nested universe is
  !! in the centre of the lattice cell). Optionally an offset map can be provided, determining
  !! whether to apply an offset in a given cell position. This can disable the local universe
  !! offset. Alternatively a single offset flag can be provided, disabling offset in all cells.
  !!
  !! Minimum pitch is set to 10 * SURF_TOL
  !!
  !! Sample Input Dictionary (3D):
  !!   hlatt { id 7;
  !!           type zHexLattice;
  !!           #origin (0.0 0.0 0.0); #
  !!           #rotation (0.0 0.0 0.0); #
  !!           orientation 1;
  !!           pitch (1.5 3.0);       // (radial pitch, axial pitch)
  !!           shape (3 3 2);         // (N1, N2, N3)
  !!           padMat void;
  !!           map ( 1 2 3            // Top layer, back row
  !!                 4 5 6            // Top layer, middle row
  !!                 7 8 9            // Top layer, front row
  !!
  !!                10 11 12
  !!                13 14 15
  !!                16 17 18  );
  !!   }
  !!
  !! Sample Input Dictionary (2D):
  !!   hlatt2D { id 8;
  !!             type zHexLattice;
  !!             orientation 2;
  !!             shape (2 2 0);       // 0 indicates infinite extent along the axis
  !!             pitch (1.0 0.0);     // Axial pitch is ignored, use 0.0 for clarity
  !!             padMat void;
  !!             map ( 1 2
  !!                   2 1);
  !!   }
  !!
  !! NOTE: Input in MAP for a single layer is WYSIWYG, as in latUniverse: the first row of
  !!   the map is the row with the highest j (2nd in-plane oblique co-ordinate).
  !!
  !! Private Members:
  !!   lat        -> Hexagonal lattice object which does the actual work
  !!   outline    -> Box type surface, bounding (from outside) the valid extent of the lattice
  !!   offset     -> Flag to disable all offsets
  !!   offsetMap  -> Map determining which cells have a lattice offset or not
  !!
  !! Interface:
  !!   universe interface
  !!
  type, public, extends(universe) :: hexUniverse
    private
    type(hexagonalLattice)                       :: lat
    type(box)                                    :: outline
    logical(defBool)                             :: offset = .true.
    integer(shortInt), dimension(:), allocatable :: offsetMap
  contains
    ! Superclass procedures
    procedure :: init
    procedure :: kill
    procedure :: findCell
    procedure :: distance
    procedure :: cross
    procedure :: cellOffset
    procedure :: getNormal
  end type hexUniverse

contains

  !!
  !! Initialise Universe
  !!
  !! See universe_inter for details.
  !!
  !! Errors:
  !!   fatalError if input is invalid.
  !!
  subroutine init(self, fill, dict, cells, surfs, mats)
    class(hexUniverse), intent(inout)                         :: self
    integer(shortInt), dimension(:), allocatable, intent(out) :: fill
    class(dictionary), intent(in)                             :: dict
    type(cellShelf), intent(inout)                            :: cells
    type(surfaceShelf), intent(inout)                         :: surfs
    type(charMap), intent(in)                                 :: mats
    integer(shortInt), dimension(:), allocatable   :: tempI
    integer(shortInt), dimension(3)                :: sizeN
    integer(shortInt)                              :: N, i, j, outFill, val
    type(dictionary)                               :: tempDict
    integer(shortInt), dimension(:,:), allocatable :: tempMap
    character(nameLen)                             :: name
    character(100), parameter :: Here = 'init (hexUniverse_class.f90)'

    ! Setup the base class
    ! With: id, origin, rotation...
    call self % setupBase(dict)

    ! Perform offsets on every cell?
    call dict % getOrDefault(self % offset, 'offset', .true.)

    ! Initialise lattice
    call self % lat % init(dict)
    sizeN = self % lat % getSize()

    ! Build outline box
    ! Unlike a Cartesian lattice, this is a safe (padded) superset of the true extent -- see
    ! hexagonalLattice % getOutlineHalfwidth for details.
    call tempDict % init(4)
    call tempDict % store('type', 'box')
    call tempDict % store('id', 1)
    call tempDict % store('origin', self % lat % getOrigin())
    call tempDict % store('halfwidth', self % lat % getOutlineHalfwidth())
    call self % outline % init(tempDict)
    call tempDict % kill()

    ! Construct fill array
    call dict % get(tempI, 'map')

    ! Ensure size matches sizeN
    if (size(tempI) /= product(sizeN)) call fatalError(Here, &
            'Lattice map size not equal to size implied by shape. Respectively: '//&
            numToChar(size(tempI))//' '//numToChar(product(sizeN)))

    ! Flip array up-down for more natural input
    ! Reshape into rank 2 array
    tempMap = reshape(tempI, [sizeN(1), sizeN(2) * sizeN(3)])
    N = size(tempMap, 2)
    do i = 1, N/2
      call swap(tempMap(:,i), tempMap(:,N - i + 1))
    end do

    ! Find background fill and change to tempMap to uniID
    tempMap = -tempMap
    call dict % get(name, 'padMat')
    outFill = charToFill(name, mats, Here)

    ! Build fill array
    allocate(fill(self % lat % getOutID()))
    N = size(tempMap, 1)
    do j = 1, size(tempMap, 2)
      do i = 1, N
        fill(i + (j-1) * N) = tempMap(i, j)
      end do
    end do
    fill(self % lat % getOutID()) = outFill
    deallocate(tempI, tempMap)

    ! Check whether there is an offset map
    if (dict % isPresent('offsetMap')) then

      if (.not. self % offset) call fatalError(Here, 'Cannot have both an offset map '//&
              'and no offset.')

      call dict % get(tempI, 'offsetMap')

      ! Ensure size matches sizeN
      if (size(tempI) /= product(sizeN)) call fatalError(Here, &
            'Offset map size not equal to size implied by shape. Respectively: '//&
            numToChar(size(tempI))//' '//numToChar(product(sizeN)))

      ! Flip array up-down for more natural input
      ! Reshape into rank 2 array
      tempMap = reshape(tempI, [sizeN(1), sizeN(2) * sizeN(3)])
      N = size(tempMap, 2)
      do i = 1, N/2
        call swap(tempMap(:,i), tempMap(:,N - i + 1))
      end do

      allocate(self % offsetMap(product(sizeN) + 1))
      N = size(tempMap, 1)
      do j = 1, size(tempMap, 2)
        do i = 1, N
          val = tempMap(i,j)
          self % offsetMap(i + (j-1) * N) = val

          ! Check that the entries are valid
          if (val /= local .and. val /= noOffset) call fatalError(Here,&
                  'Invalid entry to the offset map. Must be one of '//numToChar(local)//&
                  ' '//numToChar(noOffset)//'. Contains: '//numToChar(val))
        end do
      end do
      ! Add an entry for the padMat
      self % offsetMap(self % lat % getOutID()) = noOffset

    end if

  end subroutine init

  !!
  !! Find local cell ID given a point
  !!
  !! See universe_inter for details.
  !!
  subroutine findCell(self, localID, cellIdx, r, u)
    class(hexUniverse), intent(inout)       :: self
    integer(shortInt), intent(out)          :: localID
    integer(shortInt), intent(out)          :: cellIdx
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u

    localID = self % lat % findCell(r, u)
    cellIdx = 0

  end subroutine findCell

  !!
  !! Return distance to the next boundary between local cells in the universe
  !!
  !! See universe_inter for details.
  !!
  subroutine distance(self, d, surfIdx, coords)
    class(hexUniverse), intent(inout) :: self
    real(defReal), intent(out)        :: d
    integer(shortInt), intent(out)    :: surfIdx
    type(coord), intent(in)           :: coords

    ! Catch case if particle is outside the lattice
    if (coords % localID == self % lat % getOutID()) then

      if (self % outline % evaluate(coords % r) > ZERO) then
        ! Genuinely outside the (padded) bounding box: safe to take one large jump to its
        ! surface, since the box is guaranteed to enclose the whole lattice.
        surfIdx = OUTLINE_SURF
        d = self % outline % distance(coords % r, coords % dir)
        return
      end if

      ! Otherwise: still background, but already inside the bounding box, i.e., in the gap
      ! between the (padded) box and the true edge of the hexagonal tiling. The lattice's own
      ! distance() resolves the nearest tiling cell directly from position, so it handles
      ! this case safely. This is because hexagonalLattice assumes an infinite tiling.
    end if

    call self % lat % distance(d, surfIdx, coords % r, coords % dir)

  end subroutine distance

  !!
  !! Cross between local cells
  !!
  !! See universe_inter for details.
  !!
  subroutine cross(self, coords, surfIdx)
    class(hexUniverse), intent(inout) :: self
    type(coord), intent(inout)        :: coords
    integer(shortInt), intent(in)     :: surfIdx

    call self % findCell(coords % localID, coords % cellIdx, coords % r, coords % dir)

  end subroutine cross

  !!
  !! Return offset for the current cell
  !!
  !! See universe_inter for details.
  !!
  function cellOffset(self, coords) result (offset)
    class(hexUniverse), intent(in) :: self
    type(coord), intent(in)        :: coords
    real(defReal), dimension(3)    :: offset
    logical(defBool)               :: doOffset

    if (allocated(self % offsetMap)) then
      doOffset = self % offsetMap(coords % localID) == local
    else
      doOffset = self % offset
    end if

    if (doOffset) then
      offset = self % lat % getOffset(coords % localID)
    else
      offset = ZERO
    end if

  end function cellOffset

  !!
  !! Return normal for the given surfIdx at a given point
  !!
  !! See universe_inter for details.
  !!
  function getNormal(self, surfIdx, coords) result (normal)
    class(hexUniverse), intent(in) :: self
    integer(shortInt), intent(in)  :: surfIdx
    type(coord), intent(in)        :: coords
    real(defReal), dimension(3)    :: normal

    select case(surfIdx)
      case(OUTLINE_SURF)
        normal = self % outline % normal(coords % r, coords % dir)

      case default
        normal = self % lat % getNormal(surfIdx)

    end select

  end function getNormal

  !!
  !! Return to uninitialised state
  !!
  subroutine kill(self)
    class(hexUniverse), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Kill local
    call self % lat % kill()
    call self % outline % kill()
    if (allocated(self % offsetMap)) deallocate(self % offsetMap)
    self % offset = .true.

  end subroutine kill

end module hexUniverse_class
