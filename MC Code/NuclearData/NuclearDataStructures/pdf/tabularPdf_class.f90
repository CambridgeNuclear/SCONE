module tabularPdf_class

  use numPrecision
  use genericProcedures, only : fatalError, searchError, linearFloorIdxClosed_Real, interpolate,&
                                isSorted
  use endfConstants

  implicit none
  private

  integer(shortInt),parameter  :: histogram  = tabPdfHistogram, &
                                  linLin     = tabPdfLinLin
  real(defReal),parameter      :: tolerance = 1.0e-23

  interface linearSearch
    module procedure linearFloorIdxClosed_Real
  end interface

  !!
  !! Simple probability table for one quantity x
  !!
  type, public :: tabularPdf
    private
    real(defReal),dimension(:),allocatable       :: x
    real(defReal),dimension(:),allocatable       :: pdf
    real(defReal),dimension(:),allocatable       :: cdf
    integer(shortInt)                            :: flag            !Interpolation flag
  contains
    generic   :: init          => initPdf, initCdf
    procedure :: sample
    procedure :: bounds
    procedure :: probabilityOf
    procedure :: kill

    procedure, private :: initPdf
    procedure, private :: initCdf

  end type tabularPdf
contains

  !!
  !! Samples x given number r in <0;1>
  !! Does not check range of r
  !! Optionaly return index of a sampled bin
  !!
  function sample(self, r, res_idx) result (x)
    class(tabularPdf),intent(in)             :: self
    real(defReal),intent(in)                 :: r    ! Random Number
    integer(ShortInt), intent(out), optional :: res_idx
    real(defReal)                            :: x
    integer(shortInt)                        :: idx
    real(defReal)                            :: f, delta, ci, pi
    character(100),parameter :: Here='sample (tabularPdf_class.f90)'

    idx = linearSearch(self % cdf,r)
    call searchError(idx,Here)

    select case (self % flag)
      case (histogram)
        ci = self % cdf(idx)
        pi = self % pdf(idx)

        x = self % x(idx) + (r-ci) / pi

      case (linLin)
        ci = self % cdf(idx)
        pi = self % pdf(idx)

        f = (self % pdf(idx+1) - self % pdf(idx)) / (self % x(idx+1) - self % x(idx))

        if (f == 0) then
          x = self % x(idx) + (r-ci) / pi
        else
          delta = sqrt( pi*pi + 2*f*(r-ci))
          x = self % x(idx) + (delta - pi) / f
        end if

      case default
        call fatalError(Here,'Unknown interpolation flag')

    end select

    ! Returned the sampled index
    if(present(res_idx)) res_idx = idx

  end function sample

  !!
  !! Subroutine assigns x(1) to x_min and x(N) to x_max
  !! Returns bounds of the probability distribution
  !!
  elemental subroutine bounds(self, x_min, x_max)
    class(tabularPdf), intent(in) :: self
    real(defReal), intent(out)    :: x_min
    real(defReal), intent(out)    :: x_max

    x_min = self % x(1)
    x_max = self % x(size(self % x))

  end subroutine bounds

  !!
  !! Returns probability of x
  !!
  function probabilityOf(self,x) result (prob)
    class(tabularPdf), intent(in) :: self
    real(defReal), intent(in)     :: x
    real(defReal)                 :: prob
    integer(shortInt)             :: idx
    character(100),parameter      :: Here='probabilityOf (tabularPdf_class.f90)'

    idx = linearSearch(self % x, x)
    call searchError(idx,Here)

    select case (self % flag)
      case (histogram)
        prob = self % pdf(idx)

      case (linLin)
        prob = interpolate( self % x(idx)  ,  &
                            self % x(idx+1),  &
                            self % pdf(idx),  &
                            self % pdf(idx+1),&
                            x                 )
      case default
        call fatalError(Here,'Unknown interpolation flag')

    end select

  end function probabilityOf

  !!
  !! Release memory occupied by the table
  !!
  elemental subroutine kill(self)
    class(tabularPdf), intent(inout) :: self

    if(allocated(self % x))   deallocate(self % x)
    if(allocated(self % pdf)) deallocate(self % pdf)
    if(allocated(self % cdf)) deallocate(self % cdf)

  end subroutine kill

  !!
  !! Initialise table using PDF only
  !!
  subroutine initPdf(self,x,pdf,flag)
    class(tabularPdf), intent(inout) :: self
    real(defReal),dimension(:),intent(in)  :: x
    real(defReal),dimension(:),intent(in)  :: pdf
    integer(shortInt),intent(in)           :: flag ! Interpolation scheme flag
    integer(shortInt)                      :: i
    character(100),parameter               :: Here='init (tabularPdf_class.f90)'

    ! Check Input
    if( size(x) /= size(pdf)) call fatalError(Here,'PDF and x have diffrent size')

    if( .not.(isSorted(x)))   call fatalError(Here,'Provided x grid is not sorted not descending')
    if ( any( pdf < 0.0 ))    call fatalError(Here,'Provided PDF contains -ve values')

    ! Initialise Data
    if(allocated(self % x))   deallocate(self % x)
    if(allocated(self % pdf)) deallocate(self % pdf)
    if(allocated(self % cdf)) deallocate(self % cdf)

    self % x = x
    self % pdf = pdf
    allocate(self % cdf(size(pdf)))

    ! Calculate CDF
    select case (flag)
      case(histogram)
        self % flag = histogram
        self % cdf(1) = 0.0
        do i=2,size(x)
          self % cdf(i) = pdf(i-1) * (x(i)-x(i-1)) + self % cdf(i-1)
        end do

        if (abs(self % cdf(size(x))-1.0_defReal) > tolerance) then
          call fatalError(Here,'Calculated CDF does not integrate to 1.0 within tolerance')
        end if

      case(linLin)
        self % flag = linLin
        self % cdf(1) = 0.0
        do i=2,size(x)
          self % cdf(i) = 0.5 * (pdf(i)+pdf(i-1)) * (x(i)- x(i-1)) + self % cdf(i-1)
        end do

        if (abs(self % cdf(size(x))-1.0_defReal) > tolerance) then
          print *, self % x
          print *, self % pdf
          print *, self % cdf
          call fatalError(Here,'Calculated CDF does not integrate to 1.0 within tolerance')
        end if

      case default
        call fatalError(Here, 'Unrecognised interpolation flag')

    end select

  end subroutine initPdf

  !!
  !! Initialise table using both PDF and CDF
  !!
  subroutine initCdf(self,x,pdf,cdf,flag)
    class(tabularPdf), intent(inout) :: self
    real(defReal),dimension(:),intent(in)  :: x
    real(defReal),dimension(:),intent(in)  :: pdf
    real(defReal),dimension(:),intent(in)  :: cdf
    integer(shortInt),intent(in)           :: flag ! Interpolation scheme flag
    character(100),parameter               :: Here='init (tabularPdf_class.f90)'

    ! Check Input
    if( size(x) /= size(pdf)) call fatalError(Here,'PDF and x have diffrent size')
    if( size(x) /= size(cdf)) call fatalError(Here,'CDF and x have diffrent size')

    if( .not.(isSorted(x)))   call fatalError(Here,'Provided x grid is not sorted not decending')
    if( .not.(isSorted(cdf))) call fatalError(Here,'Provided CDF is not sorted not descending')

    if ( any( pdf < 0.0 ))    call fatalError(Here,'Provided PDF contains -ve values')
    if ( any( cdf < 0.0 ))    call fatalError(Here,'Provided CDF contains -ve values')

    if( abs(cdf(1)) > tolerance ) call fatalError(Here,'Provided CDF does not begin with 0')

    if( abs(cdf(size(cdf))-1.0_defReal) > tolerance) then
      call fatalError(Here,'Provided CDF does not end with 1')
    end if

    ! Initialise Data
    if(allocated(self % x))   deallocate(self % x)
    if(allocated(self % pdf)) deallocate(self % pdf)
    if(allocated(self % cdf)) deallocate(self % cdf)

    self % x   = x
    self % pdf = pdf
    self % cdf = cdf

    select case (flag)
      case(histogram)
        self % flag = histogram

      case(linLin)
        self % flag = linLin

      case default
        call fatalError(Here, 'Unrecognised interpolation flag')

    end select

  end subroutine initCdf

end module tabularPdf_class
