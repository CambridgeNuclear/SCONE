module tabularrelease_class

  use numPrecision

  implicit none
  private

  interface tabularrelease
    module procedure new_tabularrelease
  end interface


  type, public :: tabularrelease
    private

  contains
    procedure :: method_name      
  end type tabularrelease

contains

  function new_tabularrelease(data) result (new)
    real(defReal), intent(in)        :: data
    type(tabularrelease),pointer      :: new


  end function new_tabularrelease
    
end module tabularrelease_class
