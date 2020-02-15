module materialHandle_inter

  use numPrecision

  implicit none
  private

  !!
  !! Material Handle allows to interact with diffrent types of materials
  !!
  !! This is top abstract class for all materials.
  !! Diffrent types of materials like MG Neutron, CE Neutron etc. are its subclasses
  !!
  !! Interface:
  !!   kill -> returns to uninitialised state
  !!
  type, public, abstract :: materialHandle
  contains
    procedure(kill),deferred :: kill
  end type materialHandle

  abstract interface
    !!
    !! Return to uninitialised state
    !!
    elemental subroutine kill(self)
      import :: materialHandle
      class(materialHandle), intent(inout) :: self

    end subroutine kill
  end interface
end module materialHandle_inter
