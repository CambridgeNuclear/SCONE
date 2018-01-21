module tabularEnergy_class

  use numPrecision
  use genericProcedures,   only : fatalError
  use tabularPdf_class,    only : tabularPdf
  use RNG_class,           only : RNG

  implicit none
  private

  interface tabularEnergy
    module procedure new_tabularEnergy
  end interface

  type, public:: tabularEnergy
    private
    type(tabularPdf) :: pdf
  contains
    procedure :: sample
    procedure :: probabilityOf

    procedure :: init

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


  subroutine init(self,E,PDF,interFlag)
    class(tabularEnergy), intent(inout)   :: self
    real(defReal),dimension(:),intent(in) :: E,PDF
    integer(shortInt),intent(in)          :: interFlag
    character(100),parameter              :: Here='init (tabularEnergy_class.f90)'


    if(size(E) /= size(PDF)) call fatalError(Here,'E and PDF have diffrent size')

    call self % pdf % init(E,PDF,interFlag)

  end subroutine init


  function new_tabularEnergy(E,PDF,interFlag) result (new)
    real(defReal),dimension(:),intent(in) :: E,PDF
    integer(shortInt),intent(in)          :: interFlag
    type(tabularEnergy),pointer           :: new

    allocate(new)
    call new % init(E,PDF,interFlag)

  end function new_tabularEnergy

end module tabularEnergy_class
