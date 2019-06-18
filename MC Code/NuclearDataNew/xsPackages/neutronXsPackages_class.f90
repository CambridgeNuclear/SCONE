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
  type, public :: neutronMacroXSs
    real(defReal) :: total            = ZERO
    real(defReal) :: elasticScatter   = ZERO
    real(defReal) :: inelasticScatter = ZERO
    real(defReal) :: capture          = ZERO
    real(defReal) :: fission          = ZERO
    real(defReal) :: nuFission        = ZERO
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

end module neutronXsPackages_class
