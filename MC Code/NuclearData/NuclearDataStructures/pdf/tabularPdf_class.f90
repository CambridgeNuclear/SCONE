module tabularPdf_class

  use numPrecision
  use genericProcedures, only : fatalError, searchError, linearFloorIdxClosed_Real, interpolate,&
                                isSorted

  implicit none
  private

  integer(shortInt),parameter  :: histogram  = 0, &
                                  linLin     = 1
  real(defReal),parameter      :: tolerance = 1.0e-23

  interface linearSearch
    module procedure linearFloorIdxClosed_Real
  end interface


  type, public :: tabularPdf
    private
    real(defReal),dimension(:,:),pointer     :: data => null()
    real(defReal),dimension(:),pointer       :: x    => null()
    real(defReal),dimension(:),pointer       :: pdf  => null()
    real(defReal),dimension(:),pointer       :: cdf  => null()
    integer(shortInt)                        :: flag            !Interpolation flag
  contains
    procedure :: init
    procedure :: sample
    procedure :: probabilityOf

  end type tabularPdf
contains

  function sample(self,r) result (x)
    class(tabularPdf),intent(in) :: self
    real(defReal),intent(in)     :: r    ! Random Number
    real(defReal)                :: x
    integer(shortInt)            :: idx
    real(defReal)                :: f, delta, ci, pi
    character(100),parameter     :: Here='sample (tabularPdf_class.f90)'

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

  end function sample

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



  subroutine init(self,x,pdf,flag)
    class(tabularPdf), intent(inout) :: self
    real(defReal),dimension(:),intent(in)  :: x
    real(defReal),dimension(:),intent(in)  :: pdf
    integer(shortInt),intent(in)           :: flag ! Interpolation scheme flag
    integer(shortInt)                      :: i
    character(100),parameter               :: Here='init (tabularPdf_class.f90)'


    if( size(x) /= size(pdf)) call fatalError(Here,'PDF and x have diffrent size')

    if( .not.(isSorted(x))) call fatalError(Here,'Provided x grid is not sorted ascending')

    if (associated(self % data)) deallocate(self % data)

    allocate (self % data(size(x),3))
    self % data(:,1) = x
    self % data(:,2) = pdf
    self % data(:,3) = 0.0

    self % x   => self % data(:,1)
    self % pdf => self % data(:,2)
    self % cdf => self % data(:,3)

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
          call fatalError(Here,'Calculated CDF does not integrate to 1.0 within tolerance')
        end if

      case default
        call fatalError(Here, 'Unrecognised interpolation flag')

    end select

  end subroutine
    
end module tabularPdf_class
