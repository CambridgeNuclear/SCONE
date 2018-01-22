module tabularMiu_class

  use numPrecision
  use genericProcedures, only: fatalError, searchError, linearFloorIdxClosed_Real, interpolate
  use miuEndfPdf_class,  only: miuEndfPdf
  use RNG_class ,        only: RNG

  implicit none
  private

  integer(shortInt),parameter  :: histogram   = 0, &
                                  linLinInter = 1

  interface tabularMiu
    module procedure new_tabularMiu
  end interface

  interface linearSearch
    module procedure linearFloorIdxClosed_Real
  end interface

  type, public,extends(miuEndfPdf) :: tabularMiu
    private
    real(defReal),dimension(:),allocatable  :: miu
    real(defReal),dimension(:),allocatable  :: PDF
    real(defReal),dimension(:),allocatable  :: CDF
    integer(shortInt)                       :: interFlag
  contains
    procedure :: sample
    procedure :: probabilityOf

    procedure,private :: init

  end type tabularMiu

contains

  function sample(self,rand) result (miu)
    class(tabularMiu), intent(in)   :: self
    class(RNG), intent(inout)       :: rand
    real(defReal)                   :: miu
    integer(shortInt)               :: idx
    real(defReal)                   :: r, f, delta, ci, pi
    character(100),parameter        :: Here='sample (tabularMiu_class.f90)'

    r = rand % get()

    idx = linearSearch(self % CDF,r)
    call searchError(idx,Here)

    select case (self % interFlag)
      case (histogram)
        ci = self % CDF(idx)
        pi = self % PDF(idx)

        miu = self % miu(idx) + (r-ci) / pi

      case (linLinInter)
        ci = self % CDF(idx)
        pi = self % PDF(idx)

        f = (self % PDF(idx+1) - self % PDF(idx)) / (self % miu(idx+1) - self % miu(idx))

        if (f == 0) then
          miu = self % miu(idx) + (r-ci) / pi
        else
          delta = sqrt( pi*pi + 2*f*(r-ci))
          miu = self % miu(idx) + (delta - pi) / f
        end if

      case default
        call fatalError(Here,'Unknown interpolation flag')

    end select

  end function sample


  function probabilityOf(self,miu) result(prob)
    class(tabularMiu), intent(in)   :: self
    real(defReal), intent(in)       :: miu
    real(defReal)                   :: prob
    integer(shortInt)               :: pdfIdx
    character(100),parameter        :: Here='probabilityOf (tabularMiu_class.f90)'

    pdfIdx = linearSearch(self % miu, miu)
    call searchError(pdfIdx,Here)

    select case (self % interFlag)
      case (histogram)
        prob = self % pdf(pdfIdx)

      case (linLinInter)
        prob = interpolate( self % miu(pdfIdx)  , &
                            self % miu(pdfIdx+1), &
                            self % pdf(pdfIdx),   &
                            self % pdf(pdfIdx+1), &
                            miu                   )

      case default
        call fatalError(Here,'Unknown interpolation flag')

    end select

  end function probabilityOf


  subroutine init(self,miu,PDF,CDF,interFlag)
    class(tabularMiu), intent(inout)      :: self
    real(defReal),dimension(:),intent(in) :: miu,PDF,CDF
    integer(shortInt),intent(in)          :: interFlag

    if(allocated(self % miu)) deallocate(self % miu)
    if(allocated(self % PDF)) deallocate(self % PDF)
    if(allocated(self % CDF)) deallocate(self % CDF)

    self % interFlag = interFlag
    self % miu = miu
    self % PDF = PDF
    self % CDF = CDF

  end subroutine


  function new_tabularMiu(miu,PDF,CDF,interFlag)
    real(defReal),dimension(:),intent(in) :: miu,PDF,CDF
    integer(shortInt),intent(in)          :: interFlag
    type(tabularMiu),pointer              :: new_tabularMiu

    allocate(new_tabularMiu)
    call new_tabularMiu % init(miu,PDF,CDF,interFlag)

  end function new_tabularMiu


end module tabularMiu_class
