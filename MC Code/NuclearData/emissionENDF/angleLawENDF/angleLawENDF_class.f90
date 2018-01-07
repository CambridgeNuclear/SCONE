module angleLawENDF_class
  implicit none
  private

  type, public :: angleLawENDF_class
      private
      integer :: field_name
    contains
      procedure :: method_name
      final :: destructor
  end type angleLawENDF_class

contains

  subroutine method_name(self)
      class(angleLawENDF_class), intent(in) :: self
  end subroutine
    
end module angleLawENDF_class
