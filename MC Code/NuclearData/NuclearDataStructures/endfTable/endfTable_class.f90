module endfTable_class

  use numPrecision
  use genericProcedures, only : binarySearch, endfInterpolate, interpolate,&
                                linearCeilingIdxOpen_shortInt,  linearFloorIdxClosed_Real, &
                                searchError, fatalError, isSorted

  implicit none
  private

  interface ceilingSearch
    module procedure linearCeilingIdxOpen_shortInt
  end interface

  interface floorSearch
    module procedure linearFloorIdxClosed_Real
  end interface

  interface endfTable
    module procedure newSimple_endfTable
    module procedure newInter_endfTable
  end interface


  type, public :: endfTable
    private
    real(defReal),dimension(:),allocatable     :: x
    real(defReal),dimension(:),allocatable     :: y

    integer(shortInt)                          :: nRegions  ! Number of interpolation regions

    integer(shortInt),dimension(:),allocatable :: bounds    ! Boundaries of ENDF interpolation regions
    integer(shortInt),dimension(:),allocatable :: interENDF ! ENDF interpolation numbers

  contains
    procedure :: at
    generic   :: init => initSimple, initInter
    procedure :: kill

    procedure, private :: initSimple
    procedure, private :: initInter
  end type endfTable

contains

  subroutine initSimple(self, x, y)
    class(endfTable), intent(out)         :: self
    real(defReal),dimension(:),intent(in) :: x,y
    character(100),parameter              :: Here='initSimple (endfTable_class.f90)'

    ! Check if x and y match and if x is sorted acending array.
    if (size(x) /= size(y)) call fatalError(Here,'x and y have diffrent size!')
    if ( .not.(isSorted(x))) call fatalError(Here,'x is not sorted ascending')

    if (allocated(self % x)) deallocate(self % x)
    if (allocated(self % y)) deallocate(self % y)

    self % x = x
    self % y = y
    self % nRegions = 0

  end subroutine initSimple


  subroutine initInter(self, x, y, bounds, interENDF)
    class(endfTable), intent(out)             :: self
    real(defReal),dimension(:),intent(in)     :: x,y
    integer(shortInt),dimension(:),intent(in) :: bounds, interENDF
    character(100),parameter                  :: Here='initInter (endfTable_class.f90)'

    ! Perform Checks
    if (size(x) /= size(y))               call fatalError(Here,'x and y have diffrent size!')

    if (size(bounds) /= size (interENDF)) call fatalError(Here, 'bounds and interENDF &
                                                                & have different size')
    if ( count( bounds <= 0.0 ) > 0 )     call fatalError(Here,'bounds has -ve values')

    if (.not.isSorted(bounds))     call fatalError(Here,'bounds is not sorted')

    if ( maxval(bounds) > size(x)) call fatalError(Here,'bounds contains values larger then size(x)')

    ! Allocate Data
    if (allocated(self % x)) deallocate(self % x)
    if (allocated(self % y)) deallocate(self % y)

    if (allocated(self % bounds)) deallocate(self % bounds)
    if (allocated(self % interENDF)) deallocate(self % interENDF)

    self % x = x
    self % y = y

    self % nRegions = size(bounds)
    self % bounds = bounds
    self % interENDF = interENDF

  end subroutine initInter

  function at(self,x) result (y)
    class(endfTable),intent(in) :: self
    real(defReal), intent(in)   :: x
    real(defReal)               :: y
    integer(shortInt)           :: x_idx
    integer(shortInt)           :: bounds_idx
    real(defReal)               :: x_1, x_0, y_1, y_0
    character(100),parameter    :: Here='at (endfTable_class.f90)'

    x_idx = floorSearch(self % x, x)
    call searchError(x_idx,Here)

    x_0 = self % x(x_idx)
    x_1 = self % x(x_idx+1)

    y_0 = self % y(x_idx)
    y_1 = self % y(x_idx+1)

    select case (self % nRegions)
      case (0)
        y = interpolate(x_0, x_1, y_0, y_1, x)
      case (1)
        y = endfInterpolate(x_0, x_1, y_0, y_1, x, self%interENDF(1))
      case default
        bounds_idx = ceilingSearch(self % bounds, x_idx+1)
        call searchError(bounds_idx,Here)

        y = endfInterpolate(x_0, x_1, y_0, y_1, x, self%interENDF(bounds_idx))
    end select

  end function at

  !!
  !! Deconstruct object
  !!
  elemental subroutine kill(self)
    class(endfTable), intent(inout) :: self

    if(allocated(self % x)) deallocate(self % x)
    if(allocated(self % y)) deallocate(self % y)
    if(allocated(self % bounds)) deallocate(self % bounds)
    if(allocated(self % interENDF)) deallocate(self % interENDF)
    self % nRegions = 0

  end subroutine kill


  function newSimple_endfTable(x, y) result (new)
    real(defReal),dimension(:),intent(in) :: x,y
    type(endfTable),pointer          :: new

    allocate(new)
    call new % init(x, y)

  end function newSimple_endfTable
    

  function newInter_endfTable(x, y, bounds, interENDF) result (new)
    real(defReal),dimension(:),intent(in)     :: x, y
    integer(shortInt),dimension(:),intent(in) :: bounds, interENDF
    type(endfTable),pointer              :: new

    allocate(new)
    call new % init(x, y, bounds, interENDF)

  end function newInter_endfTable


end module endfTable_class
