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
  !!   intagral -> Integrate table from x_min to multiple values x
  !!   reloadY -> Change values on y-axis keeping x & interpolation unchanged
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
    procedure :: integral
    procedure :: reloadY
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
  !! Change values in y-axis without any othe change
  !!
  !! This subroutine is intended to eliminate memory reallocation
  !! when ENDF table is used to work with data on 2D grid. It allows
  !! to move table along 2nd dimension by reloading y-data.
  !!
  !! Args:
  !!   y [in] -> New values on y axis
  !!
  !! Error:
  !!   FatalError if y does not match current size of the table
  !!
  subroutine reloadY(self, y)
    class(endfTable), intent(inout)         :: self
    real(defReal), dimension(:), intent(in) :: y
    character(100), parameter :: Here = 'reload endfTable_class.f90'

    if (.not.allocated(self % y)) then
      call fatalError(Here, 'Cannot reload y-values on uninitialised table.')

    else if (size(y) /= size(self % y)) then
      call fatalError(Here, 'Given y-values have size: '//numToChar(size(y))//' which is diffrent &
                            &from current size of the table '//numToChar(size(self % y)))
    end if

    self % y = y

  end subroutine reloadY


  !!
  !! Integrate table from x_min to x
  !!
  !! x is given as array to allow integration of multiple points
  !! at the same time and reduce required amount of computation
  !!
  !! Args:
  !!   x [in] -> Sorted array of points in <x_min, x_max>
  !!
  !! Result:
  !!   Array of definite integral values of table from x_min to corresponding x.
  !!   For example, for x = [2] result is y = [a] a = /int_{x_min}^{2} y(t) dt
  !!   where y(t) represents tabulated values.
  !!
  !! Error:
  !!   fatalError if x is not sorted
  !!   fatalError if any value of x is outside <x_min, x_max>
  !!
  function integral(self, x) result(y)
    class(endfTable), intent(in)            :: self
    real(defReal), dimension(:), intent(in) :: x
    real(defReal), dimension(size(x))       :: y
    integer(shortInt)                       :: i, flag, reg, val
    real(defReal)                           :: x0, x1, y0, y1
    real(defReal)                           :: csum
    character(100), parameter :: Here = 'integral (endfTable_class.f90)'

    ! Preconditions
    if (.not.isSorted(x)) then
      call fatalError(Here, 'Upper bounds of integration must be given as a sorted array!')

    else if (x(size(x)) > self % x(size(self % x))) then
      call fatalError(Here, 'Values beyond upper limit of the grid were requested')

    else if (x(1) < self % x(1)) then
      call fatalError(Here, 'Values below lower limit of the grid were requested')

    end if

    ! Interpolate
    ! reg -> current interpolation region pointer
    ! val -> current place on x grid pointer
    !
    flag = linLinInterpolation
    reg =  1
    val = 1
    csum = ZERO

    do i = 2, size(self % x)
      ! Get interpolation
      if (allocated(self % bounds)) then
        if (i > self % bounds(reg)) reg = reg + 1
        flag = self % interEndf(reg)
      end if

      ! Get top and bottom of the bin
      x1 = self % x(i)
      x0 = self % x(i-1)
      y1 = self % y(i)
      y0 = self % y(i-1)

      ! Interpolate requested values
      ! Loop until valus cross into next bin. Don't forget your terminator...
      ! Yes... I've looped over few extra arrays in 1st try...
      if (val <= size(x)) then
        do while (x(val) <= self % x(i))
          y(val) = csum + endf_bin_integral(x0, x1, y0, y1, x(val), flag)
          val = val + 1
          if (val > size(x)) exit
        end do
      end if

      ! Increase csum
      ! Need to allow two neighbouring bins having same value of x
      ! In ENDF this is used to create tables with discontinuities
      if ( x1 /= x0) then
        csum = csum + endf_bin_integral(x0, x1, y0, y1, x1, flag)
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
        if (abs(f + ONE) < 1.0E-14) then
          ! Need to avoid division by almost 0
          ! assume f == -1
          int = y_0 * x_0 * log(x/x_0)

        else
          int = y_0/((f + ONE) * x_0**f) * (x**(f + ONE) - x_0**(f + ONE))

        end if

      case default
        call fatalError(Here, 'Uknown interpolation flag: '//numToChar(flag))
        int = ZERO ! Make Compiler happy

    end select

  end function endf_bin_integral

end module endfTable_class
