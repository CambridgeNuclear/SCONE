module particle_module
  implicit none
  private

  type, public :: particle
    private
    integer :: field_name
  contains
        procedure :: method_name
  end type particle

contains

  subroutine method_name(self)
    class(particle), intent(in) :: self
  end subroutine
    
end module particle_module
