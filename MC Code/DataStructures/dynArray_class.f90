module dynArray_class

  use numPrecision
  use genericProcedures, only : fatalError

  implicit none
  private

  !! Parmeters
  integer(shortInt),parameter   :: MIN_SIZE = 5               ! Initial size
  real(defReal), parameter      :: EXPAN_RATIO = 1.5_defReal  ! Expansion ratio



  !!
  !! Dynamic array for shortIntegers
  !!  Indexes are in <1;size>
  !!
  type, public :: dynIntArray
    private
    integer(shortInt), dimension(:), allocatable :: array
    integer(shortInt)                            :: mySize = 0
  contains
    ! Public interface
    generic    :: add           => add_scalar_shortInt, add_array_shortInt
    procedure  :: resize        => resize_shortInt
    procedure  :: isEmpty       => isEmpty_shortInt
    procedure  :: capacity      => capacity_shortInt
    procedure  :: getSize       => getSize_shortInt
    procedure  :: expose        => expose_shortInt
    procedure  :: get           => get_shortInt
    generic    :: assignment(=) => assign_shortInt
    procedure  :: shrink        => shrink_shortInt
    procedure  :: pop           => pop_shortInt
    procedure  :: empty         => empty_shortInt
    procedure  :: kill          => kill_shortInt

    ! Private procedures

    ! Procedures behind generic interfaces
    procedure, private  :: add_scalar_shortInt
    procedure, private  :: add_array_shortInt
    procedure, private  :: assign_shortInt

  end type dynIntArray

contains

  !!
  !! Add entry at the end of dynamic Array
  !!
  pure subroutine add_scalar_shortInt(self, entry)
    class(dynIntArray), intent(inout) :: self
    integer(shortInt), intent(in)     :: entry
    integer(shortInt)                 :: newSize

    ! Grow array if required
    newSize = self % mySize +1
    call self % resize(newSize)

    ! Store entry
    self % array(newSize) = entry
    self % mySize = newSize

  end subroutine add_scalar_shortInt

  !!
  !! Add array entry at the end of dynamic Array
  !!
  pure subroutine add_array_shortInt(self, entry)
    class(dynIntArray), intent(inout)           :: self
    integer(shortInt), dimension(:), intent(in) :: entry
    integer(shortInt)                           :: newSize

    ! Grow array if required
    newSize = self % mySize + size(entry)
    call self % resize(newSize)

    ! Store entry
    self % array(self % mySize + 1 : newSize) = entry
    self % mySize = newSize

  end subroutine add_array_shortInt

  !!
  !! Make shure that array will have memory to fit (newSize) of entries
  !! NOTE: does nothing if newSize <= current memory size
  !!
  pure subroutine resize_shortInt(self, newSize)
    class(dynIntArray), intent(inout)          :: self
    integer(shortInt), intent(in)              :: newSize
    integer(shortInt)                          :: nextSize
    integer(shortInt),dimension(:),allocatable :: tempArray

    if(.not.allocated(self % array)) then
      allocate(self % array( max(MIN_SIZE, newSize)))

    else if( newSize > size(self % array)) then
      ! Calculate next size
      nextSize = int(EXPAN_RATIO * size(self % array))
      nextSize = max(nextSize, newSize)

      ! Make temporary copy
      allocate(tempArray(nextSize))
      tempArray(1:self % mySize) = self % array(1:self % mySize)

      ! Move allocation
      call move_alloc(tempArray, self % array)

    end if
  end subroutine resize_shortInt

  !!
  !! Return .true. if array is empty or uninitialised
  !!
  pure function isEmpty_shortInt(self) result(isIt)
    class(dynIntArray), intent(in) :: self
    logical(defBool)               :: isIt

    isIt = self % mySize == 0 .or. .not.allocated(self % array)

  end function isEmpty_shortInt
    
  !!
  !! Return current memory capacity of the dynamicArray
  !!
  pure function capacity_shortInt(self) result(capacity)
    class(dynIntArray), intent(in) :: self
    integer(shortInt)              :: capacity

    if(allocated(self % array)) then
      capacity = size(self % array)

    else
      capacity = 0

    end if

  end function capacity_shortInt

  !!
  !! Get size of the dynamic array
  !!
  pure function getSize_shortInt(self) result(S)
    class(dynIntArray), intent(in) :: self
    integer(shortInt)              :: S

    S = self % mySize

  end function getSize_shortInt

  !!
  !! Return the entire dynamic array as a static array
  !!
  function expose_shortInt(self) result(res)
    class(dynIntArray), intent(in)              :: self
    integer(shortInt), dimension(self % mySize) :: res
    character(100), parameter :: Here = 'expose_shortInt (dynArray_Class.f90)'

    if (allocated(self % array)) then
      res = self % array(1: self % mySize)
    else
      call fatalError(Here,'Array is not allocated')
    end if

  end function expose_shortInt

  !!
  !! Return the entire dynamic array as a static array
  !!
  function get_shortInt(self, idx) result(res)
    class(dynIntArray), intent(in)              :: self
    integer(shortInt), intent(in)               :: idx
    integer(shortInt)                           :: res
    character(100), parameter :: Here = 'expose_shortInt (dynArray_Class.f90)'

    if (allocated(self % array)) then
      res = self % array(idx)
    else
      call fatalError(Here,'Array is not allocated')
      res = ZERO
    end if

  end function get_shortInt

  !!
  !! Assignment from a shortInt array
  !!
  subroutine assign_shortInt(LHS, RHS)
    class(dynIntArray), intent(out)            :: LHS
    integer(shortInt),dimension(:), intent(in) :: RHS

    LHS % array = RHS
    LHS % mySize = size(RHS)

  end subroutine assign_shortInt

  !!
  !! Make memory fit current size
  !!
  subroutine shrink_shortInt(self)
    class(dynIntArray), intent(inout)          :: self
    integer(shortInt),dimension(:),allocatable :: tempArray

    if(allocated( self % array)) then
      tempArray = self % array(1:self % mySize)
      call move_alloc(tempArray, self % array)
    end if

  end subroutine shrink_shortInt

  !!
  !! Pop element from the end of dynArray
  !! Return fatalError if unallocated or empty
  !!
  function pop_shortInt(self) result(res)
    class(dynIntArray), intent(inout) :: self
    integer(shortInt)                 :: res
    character(100),parameter :: Here = 'pop_shortInt (dynArray_class.f90)'

    if(self % isEmpty()) then
      call fatalError(Here,'Poping from empty array')
    end if

    res = self % array(self % mySize)
    self % mySize = self % mySize - 1

  end function pop_shortInt

  !!
  !! Remove all elements from dynamic array without deallocating memory
  !!
  pure subroutine empty_shortInt(self)
    class(dynIntArray), intent(inout) :: self

    self % mySize = 0

  end subroutine empty_shortInt

  !!
  !! Deallocate dynamic array
  !!
  pure subroutine kill_shortInt(self)
    class(dynIntArray), intent(inout) :: self

    self % mySize = 0
    if(allocated(self % array)) deallocate(self % array)

  end subroutine kill_shortInt

end module dynArray_class
