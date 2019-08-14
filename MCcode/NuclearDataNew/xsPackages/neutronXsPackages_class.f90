!!
!! This module brakes standard rules
!! It contains a library of XS Packages for Neutron particle type
!!
!!
module neutronXsPackages_class

  use numPrecision

  implicit none
  private

  !!
  !! Neutron MACROscopic Reaction XSS
  !!
  !! Public Members:
  !!   total            -> total Cross-Section [1/cm]
  !!   elasticScatter   -> sum of MT=2 elastic neutron scattering [1/cm]
  !!   inelasticScatter -> sum of all neutron producing reaction that are not elastic scattering
  !!     or fission. [1/cm]
  !!   capture          -> sum of all reactions without secendary neutrons excluding fission [1/cm]
  !!   fission          -> total Fission MT=18 Cross-section [1/cm]
  !!   nuFission        -> total average neutron production Cross-section [1/cm]
  !!
  !!  Interface:
  !!    clean -> Set all XSs to 0.0
  !!    add   -> Add a nuclide microscopic XSs to macroscopic
  !!
  type, public :: neutronMacroXSs
    real(defReal) :: total            = ZERO
    real(defReal) :: elasticScatter   = ZERO
    real(defReal) :: inelasticScatter = ZERO
    real(defReal) :: capture          = ZERO
    real(defReal) :: fission          = ZERO
    real(defReal) :: nuFission        = ZERO
  contains
    procedure :: clean => clean_neutronMacroXSs
    procedure :: add   => add_neutronMacroXSs
  end type neutronMacroXSs


  !!
  !! Neutron microscopic Reaction XSS
  !!
  !! Public Members:
  !!   total            -> total Cross-Section [barn]
  !!   elasticScatter   -> MT=2 elastic neutron scattering [barn]
  !!   inelasticScatter -> all neutron producing reaction that are not elastic scattering
  !!     or fission. [barn]
  !!   capture          -> all reactions without secendary neutrons excluding fission [barn]
  !!   fission          -> total Fission MT=18 Cross-section [barn]
  !!   nuFission        -> total average neutron production Cross-section [barn]
  !!
  type, public :: neutronMicroXSs
    real(defReal) :: total            = ZERO
    real(defReal) :: elasticScatter   = ZERO
    real(defReal) :: inelasticScatter = ZERO
    real(defReal) :: capture          = ZERO
    real(defReal) :: fission          = ZERO
    real(defReal) :: nuFission        = ZERO
  end type neutronMicroXSs

contains

  !!
  !! Clean neutron MacroXSs
  !!
  !! Sets all XSs to 0.0
  !!
  !! Args:
  !!   None
  !!
  !! Errors:
  !!   None
  !!
  elemental subroutine clean_neutronMacroXSs(self)
    class(neutronMacroXSs), intent(inout) :: self

    self % total            = ZERO
    self % elasticScatter   = ZERO
    self % inelasticScatter = ZERO
    self % capture          = ZERO
    self % fission          = ZERO
    self % nuFission        = ZERO

  end subroutine clean_neutronMacroXSs

  !!
  !! Add nuclide XSs on Macroscopic XSs
  !!
  !! Takes microscopic XSs * density and adds them to neutronMacroXSs
  !!
  !! Args:
  !!   micro [in] -> microscopic XSs
  !!   dens  [in] -> nuclide density in [1/barn/cm]
  !!
  !! Errors:
  !!   None
  !!
  elemental subroutine add_neutronMacroXSs(self, micro, dens)
    class(neutronMacroXSs), intent(inout) :: self
    type(neutronMicroXSs), intent(in)     :: micro
    real(defReal), intent(in)             :: dens

    self % total            = self % total            + dens * micro % total
    self % elasticScatter   = self % elasticScatter   + dens * micro % elasticScatter
    self % inelasticScatter = self % inelasticScatter + dens * micro % inelasticScatter
    self % capture          = self % capture          + dens * micro % capture
    self % fission          = self % fission          + dens * micro % fission
    self % nuFission        = self % nuFission        + dens * micro % nuFission

  end subroutine add_neutronMacroXSs


end module neutronXsPackages_class
