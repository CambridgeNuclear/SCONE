module tabularMu_class

  use numPrecision
  use genericProcedures, only: fatalError
  use muEndfPdf_class,   only: muEndfPdf
  use RNG_class ,        only: RNG
  use tabularPdf_class,  only: tabularPdf

  implicit none
  private



  interface tabularMu
    module procedure new_tabularMu
    module procedure new_tabularMu_withCDF
  end interface


  type, public,extends(muEndfPdf) :: tabularMu
    private
    type(tabularPdf) :: pdf
  contains
    procedure :: sample
    procedure :: probabilityOf

    generic,private   :: init => init_withPDF, init_withCDF
    procedure,private :: init_withPDF
    procedure,private :: init_withCDF


  end type tabularMu

contains

  function sample(self,rand) result (mu)
    class(tabularMu), intent(in)    :: self
    class(RNG), intent(inout)       :: rand
    real(defReal)                   :: mu
    real(defReal)                   :: r

    r = rand % get()
    mu = self % pdf % sample(r)

  end function sample


  function probabilityOf(self,mu) result(prob)
    class(tabularMu), intent(in)    :: self
    real(defReal), intent(in)       :: mu
    real(defReal)                   :: prob

    prob = self % pdf % probabilityOf(mu)

  end function probabilityOf


  subroutine init_withPDF(self,mu,PDF,interFlag)
    class(tabularMu), intent(inout)       :: self
    real(defReal),dimension(:),intent(in) :: mu,PDF
    integer(shortInt),intent(in)          :: interFlag
    character(100),parameter              :: Here='init (tabularMu_class.f90)'

    ! Check if the first element of the mu grid corresponds to -1 and last to 1
    if ( (mu(1) /= -1.0_defReal) .or. (mu(size(mu)) /= 1.0_defReal)) then

       call fatalError(Here, 'Provided mu does not begin with -1 and ends with 1')

    end if


    call self % pdf % init(mu,PDF,interFlag)

  end subroutine init_withPDF


  subroutine init_withCDF(self,mu,PDF,CDF,interFlag)
    class(tabularMu), intent(inout)       :: self
    real(defReal),dimension(:),intent(in) :: mu,PDF,CDF
    integer(shortInt),intent(in)          :: interFlag
    character(100),parameter              :: Here='init (tabularMu_class.f90)'

    ! Check if the first element of the mu grid corresponds to -1 and last to 1
    if ( (mu(1) /= -1.0_defReal) .or. (mu(size(mu)) /= 1.0_defReal)) then

       call fatalError(Here, 'Provided mu does not begin with -1 and ends with 1')

    end if


    call self % pdf % init(mu,PDF,CDF,interFlag)

  end subroutine init_withCDF


  function new_tabularMu(mu,PDF,interFlag)
    real(defReal),dimension(:),intent(in) :: mu,PDF
    integer(shortInt),intent(in)          :: interFlag
    type(tabularMu),pointer               :: new_tabularMu

    allocate(new_tabularMu)
    call new_tabularMu % init(mu,PDF,interFlag)

  end function new_tabularMu


  function new_tabularMu_withCdf(mu,PDF,CDF,interFlag) result(new)
    real(defReal),dimension(:),intent(in) :: mu,PDF,CDF
    integer(shortInt),intent(in)          :: interFlag
    type(tabularMu),pointer               :: new

    allocate(new)
    call new % init(mu,PDF,CDF,interFlag)

  end function new_tabularMu_withCDF

end module tabularMu_class
