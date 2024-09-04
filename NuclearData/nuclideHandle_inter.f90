module nuclideHandle_inter

  use numPrecision

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
  !!   getMass -> returns the mass of the nuclide (in neutron masses)
  !!   getkT -> returns nuclide temperature [MeV] of the data stored
  !!   kill -> returns to uninitialised state
  !!
  type, public, abstract :: nuclideHandle
  contains
    procedure(kill),    deferred :: kill
    procedure(getMass), deferred :: getMass
    procedure(getkT),   deferred :: getkT
  end type nuclideHandle

  abstract interface

    !!
    !! Return a mass of the nuclide
    !!
    !! Args:
    !!   None
    !!
    !! Result:
    !!   Nuclide mass in neutron masses
    !!
    !! Errors:
    !!   None
    !!
    pure function getMass(self) result(M)
      import :: nuclideHandle, defReal
      class(nuclideHandle), intent(in) :: self
      real(defReal)                    :: M
    end function getMass

    !!
    !! Return nuclide temperature
    !!
    !! Args:
    !!   None
    !!
    !! Result:
    !!   Temperature [MeV] at which the data stored in the nuclide was evaluated
    !!
    !! Errors:
    !!   None
    !!
    pure function getkT(self) result(kT)
      import :: nuclideHandle, defReal
      class(nuclideHandle), intent(in) :: self
      real(defReal)                    :: kT
    end function getkT

    !!
    !! Return to uninitialised state
    !!
    elemental subroutine kill(self)
      import :: nuclideHandle
      class(nuclideHandle), intent(inout) :: self
    end subroutine kill
  end interface



end module nuclideHandle_inter
