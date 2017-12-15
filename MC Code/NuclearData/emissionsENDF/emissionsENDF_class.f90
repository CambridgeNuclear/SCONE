module emissionsENDF_class

  use numPrecision

  implicit none
  private

  type,abstract, public :: emissionsENDF
    private
    integer(shortInt)  :: MT
  contains
    procedure(getAngleEnergy),deferred :: getAngleEnergy
    procedure(getEmission),deferred    :: getEmission
  end type emissionsENDF

abstract interface

  subroutine getAngleEnergy(self,angle,energy )
    import :: defReal, &
              emissionsENDF
    class(emissionsENDF), intent(in) :: self
    real(defReal), intent(inout) :: angle
    real(defReal), intent(inout) :: energy
  end subroutine

  subroutine getEmission(self,emission)
    import :: defReal, &
              emissionsENDF
    class(emissionsENDF), intent(in) :: self
    real(defReal), intent(inout)     :: emission
  end subroutine

end interface
contains
end module emissionsENDF_class
