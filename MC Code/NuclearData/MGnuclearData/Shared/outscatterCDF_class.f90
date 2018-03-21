module outscatterCDF_class

  use numPrecision
  use genericProcedures, only : isSorted, fatalError

  implicit none
  private

  !!
  !! Container for histogram probability of group to group scatter
  !!
  type, public :: outscatterCDF
    private
    real(defReal),dimension(:),allocatable :: cdf
  contains
    procedure :: init
    procedure :: invert
  end type outscatterCDF

contains

  !!
  !! Initialise from a single group-to-group scatter vector
  !!
  subroutine init(self,P0)
    class(outscatterCDF), intent(inout)   :: self
    real(defReal),dimension(:),intent(in) :: P0  ! P0 group to group scattering vector
    integer(shortInt)                     :: i
    character(100),parameter              :: Here ='init (outscatterCDF_class.f90)'

    ! Perform input checks
    if(any( P0 < 0.0 )) call fatalError(Here,'P0 contains -ve values')

    ! Allocate space for local cdf and calculate entries
    allocate(self % cdf(size(P0)) )
    self % cdf = P0/sum(P0)

    do i = 2,size(P0)
      self % cdf(i) = self % cdf(i-1) + self % cdf(i)

    end do

    ! Verify that initialisation was successful
    if(.not.isSorted(self % cdf)) then
      call fatalError(Here,'CDF is not sorted non-descending after init!')

    end if

  end subroutine init

  !!
  !! Invert cdf for provided number r=<0.0,1.0>
  !!
  function invert(self,r) result(G)
    class(outscatterCDF), intent(in)   :: self
    real(defReal), intent(in)          :: r
    integer(shortInt)                  :: G
    character(100), parameter          :: Here = 'invert (outscatterCDF_class.f90)'

    if ( r < 0.0) call fatalError(Here,'Number to invert CDF is -ve')

    ! Search linearly up the cdf.
    do G=1,size(self % cdf)
      if (r < self % cdf(G)) return
    end do

    if ( r > 1.0) call fatalError(Here,'Number to invert CDF is greater then 1')

    G = size(self % cdf)

  end function invert
    
end module outscatterCDF_class
