module tabularEnergy_class

  use numPrecision
  use genericProcedures,   only : fatalError
  use tabularPdf_class,    only : tabularPdf
  use RNG_class,           only : RNG

  implicit none
  private

  interface tabularEnergy
    module procedure new_tabularEnergy
    module procedure new_tabularEnergy_withCDF
  end interface

  type, public:: tabularEnergy
    private
    type(tabularPdf) :: pdf
  contains
    procedure :: sample
    procedure :: probabilityOf

    generic           :: init => init_withPDF, init_withCDF
    procedure,private :: init_withPDF
    procedure,private :: init_withCDF

  end type tabularEnergy

contains

  function sample(self,rand) result (E)
    class(tabularEnergy), intent(in) :: self
    class(RNG), intent(inout)        :: rand
    real(defReal)                    :: E
    real(defReal)                    :: r

    r = rand % get()
    E = self % pdf % sample(r)

  end function sample


  function probabilityOf(self,E) result(prob)
    class(tabularEnergy), intent(in)   :: self
    real(defReal), intent(in)          :: E
    real(defReal)                      :: prob

    prob = self % pdf % probabilityOf(E)

  end function probabilityOf


  subroutine init_withPDF(self,E,PDF,interFlag)
    class(tabularEnergy), intent(inout)   :: self
    real(defReal),dimension(:),intent(in) :: E,PDF
    integer(shortInt),intent(in)          :: interFlag
    character(100),parameter              :: Here='init (tabularEnergy_class.f90)'

    if(count( E < 0.0 ) > 0) call fatalError(Here,'E contains -ve values')

    call self % pdf % init(E,PDF,interFlag)

  end subroutine init_withPDF


  subroutine init_withCDF(self,E,PDF,CDF,interFlag)
    class(tabularEnergy), intent(inout)   :: self
    real(defReal),dimension(:),intent(in) :: E, PDF, CDF
    integer(shortInt),intent(in)          :: interFlag
    character(100),parameter              :: Here='init (tabularEnergy_class.f90)'


    if(count( E < 0.0 ) > 0) call fatalError(Here,'E contains -ve values')

    call self % pdf % init(E,PDF,CDF,interFlag)

  end subroutine init_withCDF


  function new_tabularEnergy(E,PDF,interFlag) result (new)
    real(defReal),dimension(:),intent(in) :: E,PDF
    integer(shortInt),intent(in)          :: interFlag
    type(tabularEnergy),pointer           :: new

    allocate(new)
    call new % init(E,PDF,interFlag)

  end function new_tabularEnergy


  function new_tabularEnergy_withCDF(E,PDF,CDF,interFlag) result (new)
    real(defReal),dimension(:),intent(in) :: E, PDF, CDF
    integer(shortInt),intent(in)          :: interFlag
    type(tabularEnergy),pointer           :: new

    allocate(new)
    call new % init(E,PDF,CDF,interFlag)

  end function new_tabularEnergy_withCDF


end module tabularEnergy_class
