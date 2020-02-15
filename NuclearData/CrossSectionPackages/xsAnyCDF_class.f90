module xsAnyCDF_class

  use numPrecision
  use genericProcedures, only : fatalError

  implicit none
  private

  type, public :: xsAnyCDF
    private
    real(defReal), dimension(:), allocatable :: cdf
    integer(shortInt),dimension(:),pointer   :: MTmask
  contains
    procedure :: init
    procedure :: invert

  end type xsAnyCDF

contains

  !!
  !! Build cdf for any set of reactions channels from vector of its cross-sections. Provide pointer
  !! to mask translating the index to MT number. Pointer is used to save space when the same mask
  !! is used by multiple CDFs (which is almost always going to be the case)
  !!
  subroutine init(self,xs,MTmask)
    class(xsAnyCDF), intent(inout)                    :: self
    real(defReal),dimension(:),intent(in)             :: xs
    integer(shortInt),dimension(:),pointer,intent(in) :: MTmask
    character(100),parameter                          :: Here = 'init (xsAnyCDF_class.f90)'
    integer(shortInt)                                 :: i, N

    ! Invalid input checks
    if(.not.associated(MTmask))  call fatalError(Here,'MT mask is not associated')
    if(size(MTmask) /= size(xs)) call fatalError(Here,'MT mask and xs vector have diffrent size')

    ! Connect pointer
    self % MTmask => MTmask

    ! Compute cdf
    N = size(xs)
    allocate(self % cdf(N))

    do i=1,size(xs)
      self % cdf(i) = sum(xs(1:i))
    end do

    self % cdf = self % cdf / self % cdf(N)

  end subroutine init

  !!
  !! Function to invert a CDF given a number r in <0.0;1.0>. MT is returned from MTmask pointer.
  !! Uses simple linear search for now. Faster methods may be implemented later
  !!
  function invert(self,r) result(MT)
    class(xsAnyCdf), intent(in)  :: self
    real(defReal), intent(in)    :: r
    integer(shortInt)            :: MT
    character(100),parameter     :: Here = 'invert (xsAnyCDF_class.f90)'
    integer(shortInt)            :: i


    ! Search the cdf with linear search
    do i = 1,size(self % cdf)
      if( r <= self % cdf(i)) then
        MT = self % MTmask(i)
        return
      end if
    end do

    call fatalError(Here,'Provided number to invert cdf must be > 1 ')

    ! Avoid compiler warning (Serpent territory?...)
    MT = 0

  end function invert

    
end module xsAnyCDF_class
