module cartesianField_class

  use numPrecision
  use universalVariables
  use genericProcedures,        only : fatalError, numToChar, swap
  use field_inter,              only : field
  use pieceConstantField_inter, only : pieceConstantField, kill_super => kill
  use coord_class,              only : coordList
  use particle_class,           only : particle
  use dictionary_class,         only : dictionary
  use box_class,                only : box
  use materialMenu_mod,         only : mm_matIdx => matIdx
  use cartesianLattice_class,   only : cartesianLattice

  implicit none
  private
  
  !!
  !! Public Pointer Cast
  !!
  public :: cartesianField_TptrCast

  integer(shortInt), parameter :: ALL_MATS = -1

  !!
  !! Piecewise constant field constructed from a lattice-like grid. 
  !! Values of the field are piecewise constant.
  !!
  !! Uses the Cartesian lattice to define the grid and perform most operations.
  !! Centre is placed at origin.
  !!
  !! As well as space, can include materials: values can be set on a coarse grid and differentiate
  !! between materials within a given grid cell. 
  !! If applying the values uniformly to all materials, can use the keyword 'all', 
  !! i.e., materials (all);
  !!
  !! Example dictionary:
  !!
  !! myField {
  !!   type cartesianField;
  !!   origin (x0 y0 z0);
  !!   shape (Nx Ny Nz);
  !!   pitch (Px Py Pz);
  !!   materials (fuel coolant);
  !!   ! Make up to size Nx * Ny * Nz, ascending first in x, then y, then z
  !!   fuel (
  !!    100 92 3.14 ...
  !!   ); 
  !!   coolant (
  !!    7 6 -2 ...
  !!   );
  !!   # default 8.0; #
  !!
  !! }
  !!
  type, public, extends(pieceConstantField) :: cartesianField
    private
    type(cartesianLattice) :: lattice
    type(box)              :: outline
    
    integer(shortInt)                            :: nMat = 0
    integer(shortInt), dimension(:), allocatable :: matIdxs
  contains
    ! Superclass procedures
    procedure :: init_dict
    procedure :: at
    procedure :: atP
    procedure :: distance
    procedure :: map
    procedure :: kill
  end type cartesianField

contains

  !!
  !! Initialisation
  !!
  subroutine init_dict(self, dict)
    class(cartesianField), intent(inout)          :: self
    class(dictionary), intent(in)                 :: dict
    type(dictionary)                              :: tempDict
    integer(shortInt)                             :: N, i, j, k, idx0
    real(defReal), dimension(:), allocatable      :: temp
    real(defReal), dimension(3)                   :: origin
    integer(shortInt), dimension(3)               :: sizeN
    real(defReal), dimension(:,:), allocatable    :: tempMap
    character(nameLen), dimension(:), allocatable :: mats
    character(100), parameter :: Here = 'init (cartesianField_class.f90)'
    
    ! Initialise the lattice
    call self % lattice % init(dict)
    sizeN = self % lattice % getSize()
    origin = self % lattice % getOrigin()

    ! Build outline box
    call tempDict % init(4)
    call tempDict % store('type', 'box')
    call tempDict % store('id', 1)
    call tempDict % store('origin', origin)
    call tempDict % store('halfwidth', abs(self % lattice % getCorner() - origin))
    call self % outline % init(tempDict)

    ! Construct fill array
    ! Detect how many materials are present
    self % nMat = dict % getSize('materials')

    call dict % get(mats, 'materials')
    
    if (any(mats == 'all') .and. self % nMat > 1) call fatalError(Here, 'Material "all" '//&
            'can only be used by itself and is a reserved name')

    allocate(self % matIdxs(self % nMat))
    if (dict % isPresent('all')) then
      self % matIdxs = ALL_MATS
    else
      do i = 1, self % nMat
        self % matIdxs(i) = mm_matIdx(mats(i))
      end do
    end if

    ! Size field value array
    self % N = product(sizeN * self % nMat)
    allocate(self % val(self % N + 1))

    ! Read field values for each material
    idx0 = 0
    do i = 1, size(mats)

      call dict % get(temp, mats(i))

      ! Flip array up-down for more natural input
      ! Reshape into rank 2 array
      tempMap = reshape(temp, [sizeN(1), sizeN(2) * sizeN(3)])
      N = size(tempMap, 2)
      do j = 1, N/2
        call swap(tempMap(:,j), tempMap(:,N - j + 1))
      end do

      ! Build fill array
      N = size(tempMap, 1)
      do j = 1, size(tempMap, 2)
        do k = 1, N
          self % val(idx0 + k + (j-1) * N) = tempMap(k, j)
        end do
      end do

      ! Increment starting position
      idx0 = idx0 + size(temp)

    end do

    ! Set default value when not in the field
    call dict % getOrDefault(self % val(self % N + 1), 'default', -INF)

  end subroutine init_dict

  !!
  !! Get value of the field at the co-ordinate point
  !!
  !! See pieceConstantField for details
  !!
  pure function at(self, coords) result(val)
    class(cartesianField), intent(in) :: self
    class(coordList), intent(in)      :: coords
    real(defReal)                     :: val
    integer(shortInt)                 :: localID
    
    localID = self % map(coords)
    val = self % val(localID)

  end function at
  
  !!
  !! Get index of the field at the co-ordinate point
  !!
  !! See pieceConstantField for details
  !!
  pure function map(self, coords) result(idx)
    class(cartesianField), intent(in) :: self
    class(coordList), intent(in)      :: coords
    integer(shortInt)                 :: idx
    integer(shortInt)                 :: idx0
    
    idx = self % lattice % findCell(coords % lvl(1) % r, coords % lvl(1) % dir)

    ! Outside the field
    if (idx == self % lattice % getOutID()) then
      idx = self % N + 1
      return

    end if

    ! Compare against material idx
    ! Ensure material is present
    if (any(self % matIdxs == coords % matIdx)) then
      idx0 = findloc(self % matIdxs, coords % matIdx, 1)
      idx = idx + (idx0 - 1) * product(self % lattice % getSize())
    else if (self % matIdxs(1) /= ALL_MATS) then
      idx = self % N + 1
    end if

  end function map
  
  !!
  !! Get value of the field at the particle's location
  !!
  !! See pieceConstantField for details
  !!
  function atP(self, p) result(val)
    class(cartesianField), intent(in) :: self
    class(particle), intent(in)       :: p
    real(defReal)                     :: val

    val = self % at(p % coords)

  end function atP
    
  !!
  !! Get distance to the next element of the field at the co-ordinate point and direction
  !!
  !! See pieceConstantField for details
  !!
  function distance(self, coords) result(d)
    class(cartesianField), intent(in) :: self
    class(coordList), intent(in)      :: coords
    real(defReal)                     :: d
    integer(shortInt)                 :: localID
    integer(shortInt)                 :: surfIdx

    ! Avoid compiler warnings
    surfIdx = 0

    localID = self % lattice % findCell(coords % lvl(1) % r, coords % lvl(1) % dir)
    
    ! Catch case if particle is outside the lattice
    if (localID == self % lattice % getOutID()) then
      d = self % outline % distance(coords % lvl(1) % r, coords % lvl(1) % dir)
      return

    end if

    ! Catch case if particle is in an excluded material
    if ((.not. any(coords % matIdx == self % matIdxs)) .and. self % matIdxs(1) /= ALL_MATS) then
      d = INF
      return
    end if

    call self % lattice % distance(d, surfIdx, localID, coords % lvl(1) % r, coords % lvl(1) % dir)
    
  end function distance
  
  !!
  !! Clean-up
  !!
  elemental subroutine kill(self)
    class(cartesianField), intent(inout) :: self
    
    ! Superclass
    call kill_super(self)

    call self % lattice % kill()
    call self % outline % kill()
    self % nMat = 0

  end subroutine kill
    
  !!
  !! Cast field pointer to cartesianField pointer
  !!
  !! Args:
  !!   source [in] -> source pointer of class field
  !!
  !! Result:
  !!   Null is source is not of cartesianField
  !!   Pointer to source if source is cartesianField type
  !!
  function cartesianField_TptrCast(source) result(ptr)
    class(field), pointer, intent(in) :: source
    type(cartesianField), pointer     :: ptr

    select type (source)
      type is (cartesianField)
        ptr => source

      class default
        ptr => null()
    end select

  end function cartesianField_TptrCast

end module cartesianField_class