module emissionsENDF_class
  implicit none
  private

  type,abstract public :: emissionsENDF
    private
    integer(shortInt)  :: MT
  contains
    procedure ::
  end type emissionsENDF

contains

  subroutine method_name(self)
    class(emissionsENDF), intent(in) :: self
  end subroutine
    
end module emissionsENDF_class
