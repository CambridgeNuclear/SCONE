!!
!! This module brakes standard rules
!! It contains a library of XS Packages for IMC particle type
!!
!! Pretty much a copy of neutronXsPackages_class, may contain unnecessary lines
!!
module IMCXsPackages_class

  use numPrecision
  use endfConstants

  implicit none
  private

  !!
  !! IMC MACROscopic Reaction XSS
  !!
  !! Public Members:
  !!   total            -> total opacity [1/cm]
  !!   elasticScatter   -> sum of MT=2 elastic IMC scattering [1/cm]
  !!   inelasticScatter -> sum of all IMC producing reaction that are not elastic scattering [1/cm]
  !!   capture          -> sum of all reactions without secondary photons [1/cm]
  !!   planck           -> frequency-normalised planck opacity of mateirial [1/cm]
  !!
  !!  Interface:
  !!    clean -> Set all XSs to 0.0
  !!    get   -> Return XS by MT number
  !!
  type, public :: IMCMacroXSs
    real(defReal) :: total            = ZERO
    real(defReal) :: elasticScatter   = ZERO
    real(defReal) :: inelasticScatter = ZERO
    real(defReal) :: capture          = ZERO
    real(defReal) :: planck           = ZERO
  contains
    procedure :: clean
    procedure :: get
  end type IMCMacroXSs

contains

  !!
  !! Clean IMC MacroXSs
  !!
  !! Sets all XSs to 0.0
  !!
  !! Args:
  !!   None
  !!
  !! Errors:
  !!   None
  !!
  elemental subroutine clean(self)
    class(IMCMacroXSs), intent(inout) :: self

    self % total            = ZERO
    self % elasticScatter   = ZERO
    self % inelasticScatter = ZERO
    self % capture          = ZERO
    self % planck           = ZERO

  end subroutine clean

  !!
  !! Return XSs by MT number
  !!
  !! Args:
  !!   MT [in] -> Requested MT number
  !!
  !! Result:
  !!   Value of the XS
  !!
  !! Errors:
  !!   Returns 0.0 for invalid MT
  !!
  elemental function get(self, MT) result(xs)
    class(IMCMacroXSs), intent(in) :: self
    integer(shortInt), intent(in)  :: MT
    real(defReal)                  :: xs

     select case(MT)
      case(macroTotal)
        xs = self % total

      case(macroCapture)
        xs = self % capture

      case(macroEscatter)
        xs = self % elasticScatter

      case(macroPlanck)
        xs = self % planck

      case default
        xs = ZERO

    end select

  end function get

end module IMCXsPackages_class
