module tabularMu_class

  use numPrecision
  use genericProcedures, only : fatalError
  use aceCard_class,     only : aceCard
  use muEndfPdf_inter,   only : muEndfPdf
  use RNG_class ,        only : RNG
  use tabularPdf_class,  only : tabularPdf

  implicit none
  private

  interface tabularMu
    module procedure new_tabularMu
    module procedure new_tabularMu_withCDF
    module procedure new_tabularMu_fromACE
  end interface

  !!
  !! Class that stores PDF of mu as a table
  !!
  type, public,extends(muEndfPdf) :: tabularMu
    private
    type(tabularPdf) :: pdf
  contains
    ! Superclass procedures
    procedure :: sample
    procedure :: probabilityOf
    procedure :: kill

    generic,private   :: init => init_withPDF, init_withCDF
    procedure,private :: init_withPDF
    procedure,private :: init_withCDF

  end type tabularMu

contains

  !!
  !! Sample mu given RNG
  !!
  function sample(self,rand) result (mu)
    class(tabularMu), intent(in)    :: self
    class(RNG), intent(inout)       :: rand
    real(defReal)                   :: mu
    real(defReal)                   :: r

    r = rand % get()
    mu = self % pdf % sample(r)

  end function sample

  !!
  !! Returns probability densioty of mu
  !! Does not check if mu is in <-1;1>
  !!
  function probabilityOf(self,mu) result(prob)
    class(tabularMu), intent(in)    :: self
    real(defReal), intent(in)       :: mu
    real(defReal)                   :: prob

    prob = self % pdf % probabilityOf(mu)

  end function probabilityOf

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(tabularMu), intent(inout) :: self

    call self % pdf % kill()

  end subroutine kill

  !!
  !! Initialise with PDF only
  !! Calculates CDF with PDF
  !!
  subroutine init_withPDF(self,mu,PDF,interFlag)
    class(tabularMu), intent(inout)       :: self
    real(defReal),dimension(:),intent(in) :: mu,PDF
    integer(shortInt),intent(in)          :: interFlag
    character(100),parameter              :: Here='init (tabularMu_class.f90)'

    ! Check if the first element of the mu grid corresponds to -1 and last to 1
    if ( (mu(1) /= -1.0_defReal) .or. (mu(size(mu)) /= 1.0_defReal)) then
       call fatalError(Here, 'Provided mu does not begin with -1 and ends with 1')

    end if

    ! Initialise Probability Table
    call self % pdf % init(mu,PDF,interFlag)

  end subroutine init_withPDF

  !!
  !! Initialise with PDF and CDF
  !!
  subroutine init_withCDF(self,mu,PDF,CDF,interFlag)
    class(tabularMu), intent(inout)       :: self
    real(defReal),dimension(:),intent(in) :: mu,PDF,CDF
    integer(shortInt),intent(in)          :: interFlag
    character(100),parameter              :: Here='init (tabularMu_class.f90)'

    ! Check if the first element of the mu grid corresponds to -1 and last to 1
    if ( (mu(1) /= -1.0_defReal) .or. (mu(size(mu)) /= 1.0_defReal)) then
       call fatalError(Here, 'Provided mu does not begin with -1 and ends with 1')

    end if

    ! Initialise Probability Table
    call self % pdf % init(mu,PDF,CDF,interFlag)

  end subroutine init_withCDF

  !!
  !! Construct table from PDF and interpolation flag
  !!
  function new_tabularMu(mu,PDF,interFlag)
    real(defReal),dimension(:),intent(in) :: mu,PDF
    integer(shortInt),intent(in)          :: interFlag
    type(tabularMu)                       :: new_tabularMu

    call new_tabularMu % init(mu,PDF,interFlag)

  end function new_tabularMu

  !!
  !! Construct table from PDF and CDF
  !!
  function new_tabularMu_withCdf(mu,PDF,CDF,interFlag) result(new)
    real(defReal),dimension(:),intent(in) :: mu,PDF,CDF
    integer(shortInt),intent(in)          :: interFlag
    type(tabularMu)                       :: new

    call new % init(mu,PDF,CDF,interFlag)

  end function new_tabularMu_withCDF

  !!
  !! Construct table from aceCard
  !! aceCard head needs to be set to beginning of data
  !! Will use CDF in ACE data
  !!
  function new_tabularMu_fromACE(ACE) result(new)
    type(aceCard), intent(inout)            :: ACE
    type(tabularMu)                         :: new
    integer(shortInt)                       :: inter
    integer(shortInt)                       :: N
    real(defReal),dimension(:),allocatable  :: mu
    real(defReal),dimension(:),allocatable  :: pdf
    real(defReal),dimension(:),allocatable  :: cdf
    character(100),parameter :: Here ='new_tabularMu_fromACE (tabularMu_class.f90)'

    ! Read sequence of data
    inter = ACE % readInt()        ! Read interpolation flag
    N     = ACE % readInt()        ! Read size of grid
    mu    = ACE % readRealArray(N) ! Read mu grid
    pdf   = ACE % readRealArray(N) ! Read pdf
    cdf   = ACE % readRealArray(N) ! Read cdf

    ! Initialise
    call new % init(mu,pdf,cdf,inter)

  end function new_tabularMu_fromACE

end module tabularMu_class
