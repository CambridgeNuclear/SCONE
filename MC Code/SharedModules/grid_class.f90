module grid_class

  use numPrecision
  use universalVariables
  use genericProcedures, only : fatalError, isSorted, binarySearch

  implicit none
  private

  !! Local parameters
  integer(shortInt), parameter :: UNDEF    = 0, &
                                  LIN      = 1, &
                                  LOGAR    = 2, &
                                  UNSTRUCT = 3


  type, public :: grid
   ! private
    integer(shortInt)                      :: type = UNDEF
    real(defReal),dimension(:),allocatable :: bins
    real(defReal)                          :: step

  contains
    ! Initialisation procedures
    generic   :: init => init_equalSpaced, init_unstruct
    procedure :: init_equalSpaced
    procedure :: init_unstruct

    ! Access procedures
    procedure :: bin
    procedure :: search
  end type grid

contains

  !!
  !! Initialise equaly linearly spaced grid
  !!
  subroutine init_equalSpaced(self, mini, maxi, N,type)
    class(grid), intent(inout)    :: self
    real(defReal), intent(in)     :: mini
    real(defReal), intent(in)     :: maxi
    integer(shortInt), intent(in) :: N
    character(*), intent(in)      :: type
    integer(shortInt)             :: i
    character(100), parameter :: Here = 'init_equalSpaced ( grid_class.f90)'

    ! Verify input
    if ( N < 1) call fatalError(Here,'Number of bins must be +ve')
    if ((maxi-mini)/maxi < FP_REL_TOL)  then
      call fatalError(Here,'Minimum value must be smaller then maximum above realtive FP tolerance')
    end if

    ! Allocate space for bin boundaries
    allocate(self % bins(N+1))

    select case(type)
      case('lin')
        ! Calculate bin width
        self % step = (maxi - mini)/N

        ! Assign bin boundaries
        self % bins(1) = mini
        do i= 2,N+1
          self % bins(i) = mini + (i-1) * self % step
        end do

        ! Select grid type
        self % type = LIN

      case('log')
        ! Check boundaries
        if (mini <= 0) call fatalError(Here,'For logarithmic grid minimum must be +ve')

        ! Calculate bin width
        self % step = log(maxi/mini)/N

        ! Assign bin boundaries
        self % bins(1) = mini
        do i=2,N+1
          self % bins(i) = self % bins(i-1) * exp(self % step)
        end do

        ! Select grid type
        self % type = LOGAR

      case default
        call fatalError(Here,'Grid type must be lin or log. Not: '//type)
    end select
  end subroutine init_equalSpaced

  !!
  !! Initialise unstructured grid
  !!
  subroutine init_unstruct(self,bins)
    class(grid), intent(inout)             :: self
    real(defReal),dimension(:), intent(in) :: bins
    character(100), parameter              :: Here = 'init_equalSpaced ( grid_class.f90)'

    ! Check that grid is sorted
    if( .not.isSorted(bins)) call fatalError(Here,'Provided grid is not sorted')

    ! Initialise
    self % bins = bins
    self % type = UNSTRUCT


  end subroutine init_unstruct

  !!
  !! Return value of the bins uder idx (lower bin limit)
  !! RETURNS -huge(bin) if idx is outside bounds
  !!
  elemental function bin(self,idx)
    class(grid), intent(in)       :: self
    integer(shortInt), intent(in) :: idx
    real(defReal)                 :: bin
    
    if (idx > 0 .and. idx <= size(self % bins)) then
      bin = self % bins(idx)
    else
      bin = -huge(bin)
    end if

  end function bin

  !!
  !! Returns index of the value in th grid
  !! Returns  valueOutsideArray if target is outside bounds
  !! Returns  valueOutsideArray if grid is uninitialised
  !!
  elemental function search(self, value) result(idx)
    class(grid), intent(in)  :: self
    real(defReal),intent(in) :: value
    integer(shortInt)        :: idx

    idx = 0
    select case(self % type)
      case(LIN)
        idx = floor((value-self % bins(1))/self % step)+1

      case(LOGAR)
        idx = floor(log(value/self % bins(1))/self % step) + 1

      case(UNSTRUCT)
        idx = binarySearch(self % bins, value)

    end select

    ! Check whether errors happend
    if(idx < 1 .or. idx >= size(self % bins)) idx = valueOutsideArray

  end function search

end module grid_class
