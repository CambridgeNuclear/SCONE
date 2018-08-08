module stack_class

  use numPrecision
  use genericProcedures, only : fatalError

  implicit none
  private

  real(defReal),parameter :: GROWTH_RATIO = 2.0

  !!
  !! Stack for short Integers
  !!
  type, public :: stackInt
    private
    integer(shortInt)                            :: top  = 0
    integer(shortInt), dimension(:), allocatable :: stack
  contains
    procedure :: push     => push_shortInt
    procedure :: pop      => pop_shortInt
    procedure :: isEmpty  => isEmpty_shortInt
    procedure :: size     => size_shortInt

    procedure, private :: resize => resize_shortInt
  end type stackInt

  !!
  !! Stack for characters of length nameLen
  !!
  type, public :: stackChar
    private
    integer(shortInt)                             :: top  = 0
    character(nameLen), dimension(:), allocatable :: stack
  contains
    procedure :: push     => push_char
    procedure :: pop      => pop_char
    procedure :: isEmpty  => isEmpty_char
    procedure :: size     => size_char

    procedure, private :: resize => resize_char
  end type stackChar


contains

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Integer stack procedures
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

! *** INT stack test
!  type(stackInt) :: stack
!  integer(shortInt) :: i, N, val
!
!  N = 300
!
!  do i=1,N
!    call stack % push(i)
!  end do
!
!  print *, stack % size()
!
!  do i=N,1,-1
!    call stack % pop(val)
!    print *, val
!  end do
!
!  print *, stack % isEmpty()
!  call stack % pop(val)


  !!
  !! Add element to the stack
  !!
  subroutine push_shortInt(self,val)
    class(stackInt), intent(inout) :: self
    integer(shortInt), intent(in)  :: val

    ! Resize if necessary
    call self % resize()

    ! Add new element
    self % top = self % top + 1
    self % stack(self % top) = val

  end subroutine push_shortInt

  !!
  !! Pop element from the stack
  !!
  subroutine pop_shortInt(self,val)
    class(stackInt), intent(inout) :: self
    integer(shortInt), intent(out) :: val
    character(100),parameter :: Here = 'pop_shortInt (stack_class.f90)'

    if( self % top == 0) call fatalError(Here,'Poping from empty stack')

    ! Remove element from the top
    val = self % stack( self % top)
    self % top = self % top - 1

  end subroutine pop_shortInt


  !!
  !! Returns .true. is stack is empty
  !!
  function isEmpty_shortInt(self) result(isIt)
    class(stackInt), intent(in) :: self
    logical(defBool)            :: isIt

    isIt = self % top == 0

  end function isEmpty_shortInt

  !!
  !! Returns current size of the stack population
  !!
  function size_shortInt(self) result(S)
    class(stackInt), intent(in) :: self
    integer(shortInt)           :: S

    S = self % top

  end function size_shortInt
    
  !!
  !! Increases size of the stack
  !!
  subroutine resize_shortInt(self)
    class(stackInt), intent(inout) :: self
    integer(shortInt)              :: S  ! Size
    integer(shortInt),dimension(:),allocatable :: temp

    ! Find size of the stack memory
    if ( allocated(self % stack)) then
      S = size(self % stack)
    else
      S = 0
    end if

    ! Extend storage space if needed
    if( S == self % top) then
      S = ceiling((S+1) * GROWTH_RATIO)
      allocate(temp(S))
      temp(1:self % top) = self % stack(1: self % top)
      call move_alloc(temp, self % stack)
    end if

  end subroutine resize_shortInt

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Character stack procedures
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
! ** stackChar test
!  implicit none
!  type(stackChar) :: stack
!  integer(shortInt) :: i, N
!  character(nameLen) :: val
!
!  N = 1000
!
!  do i=1,N
!    val = "My name is: "//numToChar(i)
!    call stack % push(val)
!  end do
!
!  print *, stack % size()
!
!  do i=N,1,-1
!    call stack % pop(val)
!    print *, val
!  end do
!
!  print *, stack % isEmpty()
!  call stack % pop(val)
!


  !!
  !! Add element to the stack
  !!
  subroutine push_char(self,val)
    class(stackChar), intent(inout) :: self
    character(nameLen), intent(in)  :: val

    ! Resize if necessary
    call self % resize()

    ! Add new element
    self % top = self % top + 1
    self % stack(self % top) = val

  end subroutine push_char

  !!
  !! Pop element from the stack
  !!
  subroutine pop_char(self,val)
    class(stackChar), intent(inout) :: self
    character(nameLen), intent(out) :: val
    character(100),parameter :: Here = 'pop_char (stack_class.f90)'

    if( self % top == 0) call fatalError(Here,'Poping from empty stack')

    ! Remove element from the top
    val = self % stack( self % top)
    self % top = self % top - 1

  end subroutine pop_char


  !!
  !! Returns .true. is stack is empty
  !!
  function isEmpty_char(self) result(isIt)
    class(stackChar), intent(in) :: self
    logical(defBool)            :: isIt

    isIt = self % top == 0

  end function isEmpty_char

  !!
  !! Returns current size of the stack population
  !!
  function size_char(self) result(S)
    class(stackChar), intent(in) :: self
    integer(shortInt)           :: S

    S = self % top

  end function size_char

  !!
  !! Increases size of the stack
  !!
  subroutine resize_char(self)
    class(stackChar), intent(inout) :: self
    integer(shortInt)               :: S  ! Size
    character(nameLen),dimension(:),allocatable :: temp

    ! Find size of the stack memory
    if ( allocated(self % stack)) then
      S = size(self % stack)
    else
      S = 0
    end if

    ! Extend storage space if needed
    if( S == self % top) then
      S = ceiling((S+1) * GROWTH_RATIO)
      allocate(temp(S))
      temp(1:self % top) = self % stack(1: self % top)
      call move_alloc(temp, self % stack)
    end if
  end subroutine resize_char

end module stack_class
