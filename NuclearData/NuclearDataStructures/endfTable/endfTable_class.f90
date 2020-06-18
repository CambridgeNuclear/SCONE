module endfTable_class

  use numPrecision
  use endfConstants
  use genericProcedures, only : endfInterpolate, interpolate,&
                                fatalError, isSorted, numToChar, &
                                ceilingSearch => linearCeilingIdxOpen_shortInt, &
                                floorSearch   =>  linearFloorIdxClosed_Real, &
                                binarySearch


  implicit none
  private

  !!
  !! Public functions
  !!
  public :: endf_bin_integral

  !!
  !! Data table supporting the ENDF specification
  !!
  !! The format of the table is described in section 0.5.2.1 of ENDF Manual:
  !!
  !! Trkov, A., M. Herman, and D. A. Brown. “ENDF-6 Formats Manual.”
  !! Data Formats and Procedures for the Evaluated Nuclear Data Files ENDF/B-VI and ENDF/B-VII,
  !! National Nuclear Data Center Brookhaven National Laboratory, Upton, NY, 2012, 11973–5000.
  !!
  !! Private Members:
  !!   x         -> x-grid of the table
  !!   y         -> y-grid of the table
  !!   nRegions  -> number of diffrent interpolation regions
  !!   bounds    -> Index of boundaries between diffrent interpolation regions
  !!   interENDF -> ENDF interpolation flag in each region
  !!
  !! Interface:
  !!   at   -> return y-value given x-value
  !!   init -> initialise from components
  !!   kill -> return to uninitalised state
  !!
  type, public :: endfTable
    private
    real(defReal),dimension(:),allocatable     :: x
    real(defReal),dimension(:),allocatable     :: y

    integer(shortInt)                          :: nRegions  = -1

    integer(shortInt),dimension(:),allocatable :: bounds
    integer(shortInt),dimension(:),allocatable :: interENDF

  contains
    procedure :: at
    generic   :: init => initSimple, initInter
    procedure :: kill

    procedure, private :: initSimple
    procedure, private :: initInter
  end type endfTable

contains

  !!
  !! Initialise simple table with one lin-lin interpolation region
  !!
  !! Args:
  !!   x [in] -> x-grid
  !!   y [in] -> y-grid
  !!
  !! Errors:
  !!   FatalError if x and y have difrent length.
  !!   FatalError if x is sorted and increasing
  !!
  subroutine initSimple(self, x, y)
    class(endfTable), intent(inout)       :: self
    real(defReal),dimension(:),intent(in) :: x
    real(defReal),dimension(:),intent(in) :: y
    character(100),parameter              :: Here='initSimple (endfTable_class.f90)'

    call self % kill()

    ! Check if x and y match and if x is sorted acending array.
    if (size(x) /= size(y))  call fatalError(Here,'x and y have diffrent size!')
    if ( .not.(isSorted(x))) call fatalError(Here,'x is not sorted increasing')

    ! Assign data
    self % x = x
    self % y = y
    self % nRegions = 0

  end subroutine initSimple

  !!
  !! Initialise complex table with multiple interpolation regions
  !! Args:
  !!   x [in] -> x-grid
  !!   y [in] -> y-grid
  !!   bounds [in] -> array of indexes ofinterpolation region boundaries
  !!   interENDF [in] -> array of ENDF interpolation flags
  !!
  !! Errors:
  !!   FatalError if x and y  or bounds and interENDF have difrent lengths.
  !!   FatalError if x is sorted and monotonic increasing
  !!   FatalError if bounds contains indexes that are not +ve or above size of x
  !!   FatalError if bounds is not sorted increasing
  !!   FatalError if last value of bounds is not equal to size of x-grid
  !!   FatalError if interENDF contains invalid ENDF interpolation flags
  !!
  subroutine initInter(self, x, y, bounds, interENDF)
    class(endfTable), intent(inout)             :: self
    real(defReal),dimension(:),intent(in)     :: x,y
    integer(shortInt),dimension(:),intent(in) :: bounds, interENDF
    character(100),parameter                  :: Here='initInter (endfTable_class.f90)'

    call self % kill()

    ! Perform Checks
    ! X and Y grid Error Cheks
    if (size(x) /= size(y))  call fatalError(Here,'x and y have diffrent size!')
    if ( .not.(isSorted(x))) call fatalError(Here,'x is not sorted increasing')

    ! Bounds and interENDF Error Checks
    if (size(bounds) /= size(interENDF)) call fatalError(Here, 'bounds and interENDF have different size')
    if ( any(bounds < 1) ) call fatalError(Here,'bounds has -ve values')
    if (.not.isSorted(bounds)) call fatalError(Here,'bounds is not sorted')
    if ( maxval(bounds) > size(x)) call fatalError(Here,'bounds contains values larger then size(x)')
    if (bounds(size(bounds)) /= size(x)) call fatalError(Here, 'Incomplete interpolation scheme.')

    ! Verify ENDF interpolation flags
    ! Assume that flags are a range of integers from histogramInterpolation to
    ! chargedParticleInterpolation
    if (any(interENDF < histogramInterpolation .or. interENDF > chargedParticleInterpolation)) then
      call fatalError(Here, 'InterENDF contains invalid ENDF interpolation flags!')
    end if

    ! Allocate Data
    self % x = x
    self % y = y

    self % nRegions = size(bounds)
    self % bounds = bounds
    self % interENDF = interENDF
  end subroutine initInter

  !!
  !! Return y-value form a table at x-value
  !!
  !! Args:
  !!   x [in] - x-value
  !!
  !! Result:
  !!   y-value at x-value
  !!
  !! Errors:
  !!   fatalError if x-value is outside range of x-grid
  !!
  function at(self, x) result (y)
    class(endfTable),intent(in) :: self
    real(defReal), intent(in)   :: x
    real(defReal)               :: y
    integer(shortInt)           :: x_idx
    integer(shortInt)           :: bounds_idx
    real(defReal)               :: x_1, x_0, y_1, y_0
    character(100),parameter    :: Here='at (endfTable_class.f90)'

    ! Find index
    x_idx = floorSearch(self % x, x)
    if( x_idx < 0) then
      call fatalError(Here,'Search of grid failed with error code:' // numToChar(x_idx))
    end if

    ! Find top and bottom of a bin in x and y grid
    x_0 = self % x(x_idx)
    x_1 = self % x(x_idx+1)

    y_0 = self % y(x_idx)
    y_1 = self % y(x_idx+1)

    ! Interpolate
    select case (self % nRegions)
      case (0) ! Simple int-int interpolation
        y = interpolate(x_0, x_1, y_0, y_1, x)

      case (1) ! Case for one interpolation region
        y = endfInterpolate(x_0, x_1, y_0, y_1, x, self % interENDF(1))

      case default ! Multiple interpolation regions
        bounds_idx = ceilingSearch(self % bounds, x_idx + 1)
        y = endfInterpolate(x_0, x_1, y_0, y_1, x, self % interENDF(bounds_idx))

    end select

  end function at

  !!
  !! Integrate
  !!
  !!
  !!
  function integral(self, x) result(int)
    class(endfTable), intent(in)            :: self
    real(defReal), dimension(:), intent(in) :: x
    real(defReal), dimension(size(x))       :: int
    integer(shortInt),dimension(size(x))    :: idxs
    integer(shortInt)                       :: i, N
    character(100), parameter :: Here = 'integral (endfTable_class.f90)'

    ! Preconditions
    if (.not.isSorted(x)) then
      call fatalError(Here, 'Upper bounds of integration must be given as a sorted array!')
    end if

    ! Find indexes
    N = size(x)
    do i =1, N
      idxs(i) = binarySearch(self % x, x(i))
      if (idxs(i) <= 0) then
        call fatalError(Here, 'Failed to find: '//numToChar(x(i))//' in the table grid')
      end if

    end do

  end function integral

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(endfTable), intent(inout) :: self

    if(allocated(self % x)) deallocate(self % x)
    if(allocated(self % y)) deallocate(self % y)
    if(allocated(self % bounds)) deallocate(self % bounds)
    if(allocated(self % interENDF)) deallocate(self % interENDF)
    self % nRegions = 0

  end subroutine kill

  !!
  !! Integrate ENDF bin
  !!
  !! Args:
  !!   x_0 [in] -> Low value on x-grid
  !!   x_1 [in] -> Top value on x-grid
  !!   y_0 [in] -> Low value on y-gird
  !!   y_1 [in] -> Top value on y-grid
  !!   x [in]   -> Top of integration (must be in <x_0;x_1> interval)
  !!   flag [in] -> ENDF interpolation flag
  !!
  !! Result:
  !!   Value of integral of y(x) from x_0 to x.
  !!   y(x) is linear under given interpolation (slope given by x_0, x_1, y_0, y_1)
  !!
  !! Error:
  !!   If x is outside <x_0; x_1> will extrapolate function and return value
  !!   of integral.
  !!
  function endf_bin_integral(x_0, x_1, y_0, y_1, x, flag) result(int)
    real(defReal), intent(in)     :: x_0
    real(defReal), intent(in)     :: x_1
    real(defReal), intent(in)     :: y_0
    real(defReal), intent(in)     :: y_1
    real(defReal), intent(in)     :: x
    integer(shortInt), intent(in) :: flag
    real(defReal)                 :: int
    real(defReal)                 :: f
    character(100), parameter :: Here = 'endf_bin_integral (endfTable_class.f90)'

    select case(flag)
      case (histogramInterpolation)
        int = (x - x_0) * y_0

      case (linLinInterpolation)
        f = (y_1 - y_0) / (x_1 - x_0)
        int = (x - x_0) * (y_0 + f/TWO * (x - x_0))

      case (logLinInterpolation)
        f = (y_1 - y_0) / log(x_1/x_0)
        int = (x - x_0) * y_0 + f * (x*(log(x/x_0) - ONE) + x_0 )

      case (linLogInterpolation)
        f = log(y_1/y_0) / (x_1 - x_0)
        int = y_0/f * (exp((x - x_0)*f) - ONE)

      case (loglogInterpolation)
        f = log(y_1/y_0) / log(x_1/x_0)
        if (abs(f + ONE) < 1.0E-10) then
          ! Need to avoid division by almost 0
          ! assume f == -1
          int = y_0 * x_0 * log(x/x_0)

        else
          int = y_0/((f + ONE) * x_0**f) * (x**(f + ONE) - x_0**(f + ONE))

        end if

      case default
        call fatalError(Here, 'Uknown interpolation flag: '//numToChar(flag))

    end select

  end function endf_bin_integral

end module endfTable_class
