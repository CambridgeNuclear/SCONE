module latUniverse_class

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
  use cartesianLattice_class, only : cartesianLattice, X_MIN, X_MAX, Y_MIN, Y_MAX, Z_MIN, Z_MAX

  implicit none
  private

  ! Parameters
  ! Note X/Y/Z MIN/MAX are provided by cartesianLattice_class.f90
  ! Outline must be different from those (which go from -1 to -6)
  ! Public so it can be accessed by unit tests
  integer(shortInt), public, parameter :: OUTLINE_SURF = -7

  ! Options for the offset map
  integer(shortInt), parameter :: local = 1, noOffset = 0

  !!
  !! 2D or 3D Cartesian lattice with constant pitch
  !!
  !! Universe consists of a lattice of fixed, finite size (e.g 17x17x2). Centre of the
  !! lattice is placed at the origin. An additional cell is placed beyond the lattice
  !! called background (or out) cell.
  !!
  !! Local ID is 1 in bottom X, Y & Z corner. It increases first with X then Y and lastly Z.
  !! Cells inside the lattice can only be filled with a universe (given as integer ID).
  !! Background cell can have any filling given by keyword (material or universe)
  !!
  !! Every lattice cell has an offset to its centre (so the centre of the nested universe
  !! is in the center of the lattice cell). Optionally an offset map can be provided, determining
  !! whether to apply an offset in a given cell position. This can disable the local universe
  !! offset. Alternatively a single offset flag can be provided, disabling offset in all cells.
  !!
  !! Minimum lattice pitch is set to 10 * SURF_TOL
  !!
  !! Sample Input Dictionary (3D):
  !!   latt { id 7;
  !!          type latUniverse;
  !!          #origin (0.0 0.0 0.0); #
  !!          #rotation (30.0 0.0 0.0); #
  !!          shape (3 2 2);
  !!          pitch (1.0 1.0 1.0);
  !!          padMat <u13>;
  !!          map (  1  2  3    // Top layer
  !!                 4  5  6    // Lower Y row
  !!                 7  8  9    // Bottom layer
  !!                10 11 12  );
  !!          #offsetMap ( 1 1 1
  !!          #            1 0 1
  !!          #            1 0 1
  !!          #            1 1 1 );
  !!          # offset 1;
  !!   }
  !!
  !! Sample Input Dictionary (2D):
  !!   latt2D { id 8;
  !!            shape (2 2 0);       // 0 indicates infinite extent in z-axis
  !!            pitch (1.0 1.0 0.0); // Any pitch is allowed in z, use 0.0 for clarity
  !!            padMat void;
  !!            map ( 1 2
  !!                  2 1);
  !!   }
  !!
  !! NOTE: Input in MAP for a single layer is WYSIWYG. Lower row in map is lower row in
  !!   geometry. There is no inversion like in other (e.g. Serpent) MC codes. Basically
  !!   input vector is flipped in Y and Z direction.
  !!
  !! Private Members:
  !!  lat        -> Cartesian lattice object which does the actual work
  !!  outline    -> Box type surface that is a boundary between lattice & background
  !!  offset     -> Flag to disable all offsets
  !!  offsetMap  -> Map determining which cells have a lattice offset or not
  !!
  !! Interface:
  !!   universe interface
  !!
  type, public, extends(universe) :: latUniverse
    private
    type(cartesianLattice)                       :: lat
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
  end type latUniverse

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
    class(latUniverse), intent(inout)                         :: self
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
    character(100), parameter :: Here = 'init (latUniverse_class.f90)'

    ! Setup the base class
    ! With: id, origin rotations...
    call self % setupBase(dict)
    
    ! Perform offsets on every cell?
    call dict % getOrDefault(self % offset, 'offset', .true.)

    ! Initialise lattice
    call self % lat % init(dict)
    sizeN = self % lat % getSize()

    ! Build outline box
    call tempDict % init(4)
    call tempDict % store('type', 'box')
    call tempDict % store('id', 1)
    call tempDict % store('origin', self % lat % getOrigin())
    call tempDict % store('halfwidth', abs(self % lat % getCorner() - self % lat % getOrigin()))
    call self % outline % init(tempDict)

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
    class(latUniverse), intent(inout)       :: self
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
    class(latUniverse), intent(inout) :: self
    real(defReal), intent(out)        :: d
    integer(shortInt), intent(out)    :: surfIdx
    type(coord), intent(in)           :: coords

    ! Catch case if particle is outside the lattice
    if (coords % localID == self % lat % getOutID()) then
      surfIdx = OUTLINE_SURF
      d = self % outline % distance(coords % r, coords % dir)
      return

    end if

    call self % lat % distance(d, surfIdx, coords % localID, coords % r, coords % dir)

  end subroutine distance

  !!
  !! Cross between local cells
  !!
  !! See universe_inter for details.
  !!
  subroutine cross(self, coords, surfIdx)
    class(latUniverse), intent(inout) :: self
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
    class(latUniverse), intent(in) :: self
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
    class(latUniverse), intent(in) :: self
    integer(shortInt), intent(in)  :: surfIdx
    type(coord), intent(in)        :: coords
    real(defReal), dimension(3)    :: normal
    character(100), parameter      :: Here = 'getNormal (latUniverse_class.f90)'

    ! Convert surfIdx to appropriate axis
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
    class(latUniverse), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Kill local
    call self % lat % kill()
    call self % outline % kill()
    if (allocated(self % offsetMap)) deallocate(self % offsetMap)
    self % offset = .true.

  end subroutine kill

end module latUniverse_class
