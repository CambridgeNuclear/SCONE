module testNeutronMaterial_class

  use numPrecision
  use particle_class,           only : particle
  use neutronMaterial_inter,    only : neutronMaterial
  use neutronXsPackages_class,  only : neutronMacroXSs

  implicit none
  private

  !!
  !! Test Neutron Material
  !! Mimics Neutron Material for testing
  !!
  !! Public Mambers
  !!   xss -> Stores Macroscopic XSs
  !!
  !! Interface:
  !!   materialHandle interface
  !!   neutronMaterial interface
  !!
  type, public, extends(neutronMaterial) :: testNeutronMaterial
    type(neutronMacroXSs) :: xss

  contains
    ! Superclass Procedures
    procedure :: kill
    procedure :: isFissile
    procedure :: getMacroXSs_byP
    procedure :: getMTxs

  end type testNeutronMaterial

contains

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(testNeutronMaterial), intent(inout) :: self

    call self % xss % clean()

  end subroutine kill

  !!
  !! Return True if material is fissile
  !!
  !! See neutronMaterial_inter for more details
  !!
  elemental function isFissile(self) result(isIt)
    class(testNeutronMaterial), intent(in) :: self
    logical(defBool)                       :: isIt

    isIt = self % xss % fission > ZERO

  end function isFissile

  !!
  !! Return Macroscopic Material XSs given particle
  !!
  !! See neutronMaterial_inter for more details
  !!
  subroutine getMacroXSs_byP(self, xss, p)
   class(testNeutronMaterial), intent(in) :: self
   type(neutronMacroXSs), intent(out)     :: xss
   class(particle), intent(in)            :: p

   xss = self % xss

  end subroutine getMacroXSs_byP

  !!
  !! Return Macroscopic XS for the material given particle and an MT number
  !!
  !! See neutronMaterial_inter for more details
  !!
  function getMTxs(self, MT, p) result(xs)
   class(testNeutronMaterial), intent(in) :: self
   integer(shortInt), intent(in)          :: MT
   class(particle), intent(in)            :: p
   real(defReal)                          :: xs

   xs = ZERO

 end function getMTxs

end module testNeutronMaterial_class
