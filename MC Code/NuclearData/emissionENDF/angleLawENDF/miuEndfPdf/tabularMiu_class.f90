module tabularMiu_class

  use numPrecision
  use miuEndfPdf_class,  only: miuEndfPdf
  use RNG_class ,        only: RNG
  use tabularPdf_class,  only: tabularPdf

  implicit none
  private



  interface tabularMiu
    module procedure new_tabularMiu
  end interface


  type, public,extends(miuEndfPdf) :: tabularMiu
    private
    type(tabularPdf) :: pdf
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
    real(defReal)                   :: r

    r = rand % get()
    miu = self % pdf % sample(r)

  end function sample


  function probabilityOf(self,miu) result(prob)
    class(tabularMiu), intent(in)   :: self
    real(defReal), intent(in)       :: miu
    real(defReal)                   :: prob

    prob = self % pdf % probabilityOf(miu)

  end function probabilityOf


  subroutine init(self,miu,PDF,interFlag)
    class(tabularMiu), intent(inout)      :: self
    real(defReal),dimension(:),intent(in) :: miu,PDF
    integer(shortInt),intent(in)          :: interFlag

    call self % pdf % init(miu,PDF,interFlag)

  end subroutine


  function new_tabularMiu(miu,PDF,interFlag)
    real(defReal),dimension(:),intent(in) :: miu,PDF
    integer(shortInt),intent(in)          :: interFlag
    type(tabularMiu),pointer              :: new_tabularMiu

    allocate(new_tabularMiu)
    call new_tabularMiu % init(miu,PDF,interFlag)

  end function new_tabularMiu


end module tabularMiu_class
