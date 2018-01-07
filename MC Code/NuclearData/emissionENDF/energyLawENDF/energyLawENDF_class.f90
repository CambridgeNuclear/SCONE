module energyLawENDF_class
  implicit none
  private

  type, public :: energyLawENDF_class
      private
      integer :: field_name
    contains
      procedure :: method_name
      final :: destructor
  end type energyLawENDF_class

contains

  subroutine method_name(self)
      class(energyLawENDF_class), intent(in) :: self
  end subroutine
    
end module energyLawENDF_class
