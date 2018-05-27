module releaseLawENDFfactory_func

  use numPrecision
  use endfConstants
  use genericProcedures,       only : fatalError
  use aceCard_class,           only : aceCard
  use releaseLawENDF_inter,    only : releaseLawENDF
  use constantRelease_class,   only : constantRelease
  use polynomialRelease_class, only : polynomialRelease
  use tabularRelease_class,    only : tabularRelease

  implicit none
  private


  integer(shortInt),parameter :: POLYNOMIAL_NU = 1 ,&
                                 TABULAR_NU    = 2


  public :: new_releaseLawENDF
  public :: new_releaseLawENDF_ptr
  public :: new_totalNU
  public :: new_promptNU
  public :: new_delayedNU

contains

  !!
  !! Returns an allocatable releaseLawENDF from aceCard and MT number
  !! aceCard can be in any poistion. Its position changes at output.
  !!
  function new_releaseLawENDF(ACE,MT) result(new)
    type(aceCard), intent(inout)      :: ACE
    integer(shortInt), intent(in)     :: MT
    class(releaseLawENDF),allocatable :: new
    integer(shortInt)                 :: TY
    character(100), parameter :: Here = 'new_releaseLawENDF ( releaseLawENDFfatory_func.f90)'

    ! Read neutron Release (TY value)
    TY = ACE % neutronReleaseMT(MT)

    select case(TY)
      case(101:) ! Special energy dependent yield -> not supported
        call fatalError(Here,'TY > 100. Spetial energy dependent yield. Not supported yet')

      case(19) ! Fission
        allocate(new, source = new_totalNu(ACE) )

      case default ! Constant release. Need to convert integer to real
        allocate(new, source = constantRelease(real(TY,defReal)) )

    end select
  end function new_releaseLawENDF


  !!
  !! Returns a pointer to allocated releaseLawENDF from aceCard and MT number
  !!
  function new_releaseLawENDF_ptr(ACE,MT) result(new)
    type(aceCard), intent(inout)  :: ACE
    integer(shortInt),intent(in)  :: MT
    class(releaseLawENDF),pointer :: new

    ! Allocate pointer and copy data from local allocatable
    allocate(new, source = new_releaseLawENDF(ACE,MT))

  end function new_releaseLawENDF_ptr

  !!
  !! Returns allocated release law for total NU data
  !! If NU data does not exist returns errror
  !!
  function new_totalNu(ACE) result(new)
    type(aceCard), intent(inout)       :: ACE
    class(releaseLawENDF), allocatable :: new
    character(100), parameter :: Here = 'new_totalNu ( releaseLawENDFfatory_func.f90)'

    ! Check if the data exists
    if(.not.ACE % hasNuTotal()) then
      call fatalError(Here, 'Total NU data is not present')

    end if

    ! Set to beginning of NU total data and allocate NU data
    call ACE % setToNuTotal()
    call allocateNu(new,ACE)

  end function new_totalNu

  !!
  !! Returns allocated release law for prompt NU data
  !! If NU data does not exist returns errror
  !! If only single prompt/total NU is provided in ACE data it is assumed to be total
  !!
  function new_promptNu(ACE) result(new)
    type(aceCard), intent(inout)       :: ACE
    class(releaseLawENDF), allocatable :: new
    character(100), parameter :: Here = 'new_promptNu ( releaseLawENDFfatory_func.f90)'

    ! Check if the data exists
    if(.not.ACE % hasNuPrompt()) then
      call fatalError(Here, 'Prompt NU data is not present')

    end if

    ! Set to beginning of NU prompt data and allocate NU data
    call ACE % setToNuPrompt()
    call allocateNu(new,ACE)

  end function new_promptNu

  !!
  !! Returns allocated release law for delayed NU data
  !! If NU data does not exist returns errror
  !!
  function new_delayedNu(ACE) result(new)
    type(aceCard), intent(inout)       :: ACE
    class(releaseLawENDF), allocatable :: new
    character(100), parameter :: Here = 'new_delayedNu ( releaseLawENDFfatory_func.f90)'

    ! Check if the data exists
    if(.not.ACE % hasNuDelayed()) then
      call fatalError(Here, 'Delayed NU data is not present')

    end if

    ! Set to beginning of NU deleyed data and allocate NU data
    call ACE % setToNuDelayed()
    call allocateNu(new,ACE)

  end function new_delayedNu

  !!
  !! Private subroutine to avoid code repeat
  !! accCard read head needs to be at the beggining of NU data block
  !!
  !!
  subroutine allocateNu(new_nu,ACE)
    class(releaseLawENDF), allocatable, intent(inout) :: new_nu
    type(aceCard), intent(inout)                      :: ACE
    integer(shortInt)                                 :: LNU
    character(100), parameter :: Here = 'allocateNu ( releaseLawENDFfatory_func.f90)'

    ! Deallocate new_nu if allocated
    if(allocated (new_nu)) deallocate(new_nu)

    ! Read type of Nu data
    LNU = ACE % readInt()

    select case(LNU)
      case(POLYNOMIAL_NU)
        allocate(new_nu, source = polynomialRelease(ACE))

      case(TABULAR_NU)
        allocate(new_nu, source = tabularRelease(ACE))

      case default
        call fatalError(Here,'Unrecoginised LNU. Not 1 or 2')

    end select
  end subroutine allocateNu


end module releaseLawENDFfactory_func
