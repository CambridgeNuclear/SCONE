module nuclideHandle_inter


  implicit none
  private

  !!
  !! Nuclide Handle allows to interact with diffrent types of nuclides
  !!
  !! This is top abstract class for all mnuclides
  !! Diffrent types of nuclides like diffrent implementations of CE Neutron data
  !! are its subclasses.
  !!
  !! Interface:
  !!   kill -> returns to uninitialised state
  !!
  type, public, abstract :: nuclideHandle
  contains
    procedure(kill), deferred :: kill
  end type nuclideHandle

  abstract interface
    !!
    !! Return to uninitialised state
    !!
    elemental subroutine kill(self)
      import :: nuclideHandle
      class(nuclideHandle), intent(inout) :: self
    end subroutine kill
  end interface


    
end module nuclideHandle_inter
