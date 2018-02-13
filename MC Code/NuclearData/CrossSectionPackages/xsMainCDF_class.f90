module xsMainCDF_class

  use numPrecision
  use endfConstants
  use genericProcedures, only : fatalError
  use xsCDF_class,       only : xsCDF

  implicit none
  private

  integer(shortInt), dimension(3), parameter :: reactionMask = [ anyScatter, anyCapture, anyFission]


  type, public,extends(xsCDF) :: xsMainCDF
   ! private
    real(defReal), dimension(3) :: cdf  ! 1- scattering 2- capture 3- fission
  contains
    procedure :: invert
    procedure :: init
    procedure :: interpolate
  end type xsMainCDF

contains

  !!
  !! Build the CDF for sampling main reaction channel from individual cross-sections
  !!
  elemental subroutine init(self,scatter,capture,fission)
    class(xsMainCDF), intent(inout)   :: self
    real(defReal),intent(in)          :: scatter
    real(defReal),intent(in)          :: capture
    real(defReal),intent(in),optional :: fission
    real(defReal)                     :: locFission
    real(defReal)                     :: total
    ! Check if the fission XS is supplied
    if(present(fission)) then
      locFission = fission
    else
      locFission = 0.0
    end if

    ! Compute CDF
    self % cdf(1) = scatter
    self % cdf(2) = scatter + capture
    self % cdf(3) = scatter + capture + locFission
    self % cdf    = self % cdf/ self % cdf(3)

  end subroutine init


  !!
  !! Invert the CDF by providing a number r in <0.0;1.0>
  !!
  function invert(self,r) result(mask)
    class(xsMainCDF), intent(in) :: self
    real(defReal), intent(in)    :: r
    integer(shortInt)            :: mask
    integer(shortInt)            :: i
    character(100),parameter     :: Here ='invert (xsMainCDF_class.f90)'

    ! Search the cdf
    do i = 1,3
      if( r <= self % cdf(i)) then
        mask = reactionMask(i)
        return
      end if
    end do

    call fatalError(Here,'Provided number to invert cdf must be > 1 ')

  end function invert

  !!
  !! Interpolate the CDF between low and top xsMainCDFs using the provided interpolation factor
  !! (x-x_0)/(x_1-x_0), where 1-top value; 0-low value. No error checking is made on output so
  !! care must be taken not to mix top with low and provide correct interpolation factor
  !!
  subroutine interpolate(self,low,top,f)
    class(xsMainCDF), intent(inout) :: self
    type(xsMainCDF), intent(in)     :: low
    type(xsMainCDF), intent(in)     :: top
    real(defReal), intent(in)       :: f
    real(defReal)                   :: f2

    ! Calculate (1-f)
    f2 = 1.0 - f
    ! Interpolate CDF
    self % cdf = f2*low % cdf + f*top % cdf

  end subroutine interpolate
    
end module xsMainCDF_class
