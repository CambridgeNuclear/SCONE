module correlatedLawENDF_class
  implicit none
  private

  type, public :: correlatedLawENDF_class
      private
      integer :: field_name
    contains
      procedure :: method_name
      final :: destructor
  end type correlatedLawENDF_class

contains

  subroutine method_name(self)
      class(correlatedLawENDF_class), intent(in) :: self
  end subroutine
    
end module correlatedLawENDF_class
