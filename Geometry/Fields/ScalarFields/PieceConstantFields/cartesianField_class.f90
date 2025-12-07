module cartesianField_class

  use numPrecision
  use universalVariables
  use genericProcedures,        only : fatalError, numToChar, swap
  use field_inter,              only : field
  use pieceConstantField_inter, only : pieceConstantField
  use coord_class,              only : coordList
  use particle_class,           only : particle
  use dictionary_class,         only : dictionary
  use box_class,                only : box
  use materialMenu_mod,         only : mm_matIdx => matIdx

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
  !! Similar to a Cartesian lattice. Centre is placed at origin.
  !! Can include materials. If applying the values uniformly to all materials, can use
  !! the keyword 'all', i.e., materials (all);
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
    real(defReal), dimension(3)     :: pitch = ZERO
    integer(shortInt), dimension(3) :: sizeN = 0
    real(defReal), dimension(3)     :: corner = ZERO
    real(defReal), dimension(3)     :: a_bar  = ZERO
    type(box)                       :: outline
    integer(shortInt)               :: outLocalID = 0
    
    integer(shortInt)                            :: nMat = 0
    integer(shortInt), dimension(:), allocatable :: matIdxs
  contains
    ! Superclass procedures
    procedure :: init
    procedure :: at
    procedure :: atP
    procedure :: distance
    procedure :: kill

    procedure, private :: getLocalID
  end type cartesianField

contains

  !!
  !! Initialisation
  !!
  subroutine init(self, dict)
    class(cartesianField), intent(inout)          :: self
    class(dictionary), intent(in)                 :: dict
    type(dictionary)                              :: tempDict
    integer(shortInt)                             :: N, i, j, k, idx0
    integer(shortInt), dimension(:), allocatable  :: tempI
    real(defReal), dimension(:), allocatable      :: temp
    real(defReal), dimension(3)                   :: origin
    real(defReal), dimension(:,:), allocatable    :: tempMap
    character(nameLen), dimension(:), allocatable :: mats
    character(100), parameter :: Here = 'init (cartesianField_class.f90)'
    
    ! Load pitch
    call dict % get(temp, 'pitch')
    N = size(temp)

    if (N /= 3) then
      call fatalError(Here, 'Pitch must have size 3. Has: '//numToChar(N))
    end if
    self % pitch = temp
    
    ! Load origin
    call dict % get(temp, 'origin')
    N = size(temp)

    if (N /= 3) then
      call fatalError(Here, 'Origin must have size 3. Has: '//numToChar(N))
    end if
    origin = temp

    ! Load Size
    call dict % get(tempI, 'shape')
    N = size(tempI)

    if (N /= 3) then
      call fatalError(Here, 'Shape must have size 3. Has: '//numToChar(N))
    else if (any(tempI < 0)) then
      call fatalError(Here, 'Shape contains -ve entries')
    end if
    self % sizeN = tempI

    ! Detect reduced Z dimension
    if (self % sizeN(3) == 0) then
      self % sizeN(3) = 1
      self % pitch(3) = TWO * INF
    end if

    ! Check X & Y for 0 size
    if (any( self % sizeN == 0)) call fatalError(Here, 'Shape in X and Y axis cannot be 0.')

    ! Check for invalid pitch
    if (any(self % pitch < 10 * SURF_TOL)) then
     call fatalError(Here, 'Pitch size must be larger than: '//numToChar( 10 * SURF_TOL))
   end if

    ! Calculate halfwidth and corner
    self % a_bar = self % pitch * HALF - SURF_TOL
    self % corner = origin -(self % sizeN * HALF * self % pitch)

    ! Build outline box
    call tempDict % init(4)
    call tempDict % store('type', 'box')
    call tempDict % store('id', 1)
    call tempDict % store('origin', origin)
    call tempDict % store('halfwidth', abs(self % corner - origin))
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
    self % outLocalID = product(self % sizeN) + 1
    self % N = product(self % sizeN * self % nMat) + 1
    allocate(self % val(self % N))

    ! Read field values for each material
    idx0 = 0
    do i = 1, size(mats)

      call dict % get(temp, mats(i))

      ! Flip array up-down for more natural input
      ! Reshape into rank 2 array
      tempMap = reshape(temp, [self % sizeN(1), self % sizeN(2) * self % sizeN(3)])
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
    call dict % getOrDefault(self % val(self % N), 'default', -INF)

  end subroutine init

  !!
  !! Get value of the field at the co-ordinate point
  !!
  !! See pieceConstantField for details
  !!
  function at(self, coords) result(val)
    class(cartesianField), intent(in) :: self
    class(coordList), intent(in)      :: coords
    real(defReal)                     :: val
    integer(shortInt)                 :: idx, localID
    
    localID = self % getLocalID(coords % lvl(1) % r, coords % lvl(1) % dir)

    if (localID == self % outLocalID) then
      val = self % val(self % N)
      return

    end if

    ! Compare against material idx
    if (self % matIdxs(1) /= ALL_MATS) then
      idx = findloc(self % matIdxs, coords % matIdx, 1)
      localID = localID + (idx - 1) * product(self % sizeN)
    end if
      
    val = self % val(localID)

  end function at
  
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
    real(defReal)                     :: test_d
    integer(shortInt)                 :: localID, temp, base, ax, i
    integer(shortInt), dimension(3)   :: ijk
    real(defReal), dimension(3)       :: bounds, r_bar, u

    localID = self % getLocalID(coords % lvl(1) % r, coords % lvl(1) % dir)
    
    ! Catch case if particle is outside the lattice
    if (localID == self % outLocalID) then
      d = self % outline % distance(coords % lvl(1) % r, coords % lvl(1) % dir)
      return

    end if
    
    ! Compute ijk of localID
    temp = localID - 1

    base = temp / self % sizeN(1)
    ijk(1) = temp - self % sizeN(1) * base + 1

    temp = base
    base = temp / self % sizeN(2)
    ijk(2) = temp - self % sizeN(2) * base + 1

    ijk(3) = base + 1

    ! Find position wrt lattice cell centre
    ! Need to use localID to properly handle under and overshoots
    u = coords % lvl(1) % dir
    r_bar = coords % lvl(1) % r - self % corner
    r_bar = r_bar - (ijk - HALF) * self % pitch

    ! Select surfaces in the direction of the particle
    bounds = sign(self % pitch * HALF, u)

    ! Find minimum distance
    ! Relay on IEEE 754 standard (for floating point numbers)
    ! 0.0/0.0 = NaN and (NaN < A = false; for every A)
    ! A/0.0 = Infinity (if A > 0.0)
    !
    ! Provide default axis to ensure no out of bounds array access if
    ! all distances happen to be infinite
    d = INF
    ax = 1
    do i = 1, 3
      ! Nominator and denominator will have the same sign (by earlier bounds selection)
      test_d = (bounds(i) - r_bar(i)) / u(i)

      if (test_d < d) then
        d = test_d
        ax = i
      end if
    end do

    ! Cap distance value
    d = max(ZERO, d)
    d = min(INF, d)

  end function distance
  
  !!
  !! Clean-up
  !!
  elemental subroutine kill(self)
    class(cartesianField), intent(inout) :: self
    
    call self % killSuper()

    self % pitch = ZERO
    self % sizeN = 0
    self % nMat = 0
    self % corner = ZERO
    self % a_bar  = ZERO
    self % outLocalID = 0
    call self % outline % kill()

  end subroutine kill

  !!
  !! Find the local integer ID in the field given position and direction
  !!
  function getLocalID(self, r, u) result(localID)
    class(cartesianField), intent(in) :: self
    real(defReal), dimension(3)       :: r
    real(defReal), dimension(3)       :: u
    integer(shortInt)                 :: localID
    real(defReal), dimension(3)       :: r_bar
    integer(shortInt), dimension(3)   :: ijk
    integer(shortInt)                 :: i, inc
          
    ijk = floor((r - self % corner) / self % pitch) + 1

    ! Get position wrt middle of the lattice cell
    r_bar = r - self % corner - ijk * self % pitch + HALF * self % pitch

    ! Check if position is within surface tolerance
    ! If it is, push it to next cell
    do i = 1, 3
      if (abs(r_bar(i)) > self % a_bar(i) .and. r_bar(i)*u(i) > ZERO) then

        ! Select increment. Ternary expression
        if (u(i) < ZERO) then
          inc = -1
        else
          inc = 1
        end if

        ijk(i) = ijk(i) + inc
      end if
    end do

    if (any(ijk <= 0 .or. ijk > self % sizeN)) then ! Point is outside lattice
      localID = self % outLocalID

    else
      localID = ijk(1) + self % sizeN(1) * (ijk(2)-1 + self % sizeN(2) * (ijk(3)-1))

    end if


  end function getLocalID
    
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
  pure function cartesianField_TptrCast(source) result(ptr)
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
