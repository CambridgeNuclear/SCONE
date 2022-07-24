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
  !!   total            -> total Cross-Section [1/cm]
  !!   elasticScatter   -> sum of MT=2 elastic IMC scattering [1/cm]
  !!   inelasticScatter -> sum of all IMC producing reaction that are not elastic scattering [1/cm]
  !!   capture          -> sum of all reactions without secondary photons [1/cm]
  !!
  !!  Interface:
  !!    clean -> Set all XSs to 0.0
  !!    add   -> Add a nuclide microscopic XSs to macroscopic
  !!    get   -> Return XS by MT number
  !!
  type, public :: IMCMacroXSs
    real(defReal) :: total            = ZERO
    real(defReal) :: elasticScatter   = ZERO
    real(defReal) :: inelasticScatter = ZERO
    real(defReal) :: capture          = ZERO
  contains
    procedure :: clean => clean_IMCMacroXSs
    procedure :: add   => add_IMCMacroXSs
    procedure :: get
    procedure :: invert => invert_macroXSs
  end type IMCMacroXSs


  !!
  !! IMC microscopic Reaction XSS
  !!
  !! Public Members:
  !!   total            -> total Cross-Section [barn]
  !!   elasticScatter   -> MT=2 elastic IMC scattering [barn]
  !!   inelasticScatter -> all photon producing reaction that are not elastic scattering [barn]
  !!   capture          -> all reactions without secendary photons [barn]
  !!
  type, public :: IMCMicroXSs
    real(defReal) :: total            = ZERO
    real(defReal) :: elasticScatter   = ZERO
    real(defReal) :: inelasticScatter = ZERO
    real(defReal) :: capture          = ZERO
  contains
    procedure :: invert => invert_microXSs
  end type IMCMicroXSs

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
  elemental subroutine clean_IMCMacroXSs(self)
    class(IMCMacroXSs), intent(inout) :: self

    self % total            = ZERO
    self % elasticScatter   = ZERO
    self % inelasticScatter = ZERO
    self % capture          = ZERO

  end subroutine clean_IMCMacroXSs

  !!
  !! Add nuclide XSs on Macroscopic XSs
  !!
  !! Takes microscopic XSs * density and adds them to IMCMacroXSs
  !!
  !! Args:
  !!   micro [in] -> microscopic XSs
  !!   dens  [in] -> nuclide density in [1/barn/cm]
  !!
  !! Errors:
  !!   None
  !!
  elemental subroutine add_IMCMacroXSs(self, micro, dens)
    class(IMCMacroXSs), intent(inout) :: self
    type(IMCMicroXSs), intent(in)     :: micro
    real(defReal), intent(in)             :: dens

    self % total            = self % total            + dens * micro % total
    self % elasticScatter   = self % elasticScatter   + dens * micro % elasticScatter
    self % inelasticScatter = self % inelasticScatter + dens * micro % inelasticScatter
    self % capture          = self % capture          + dens * micro % capture

  end subroutine add_IMCMacroXSs

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
    integer(shortInt), intent(in)      :: MT
    real(defReal)                      :: xs

     select case(MT)
      case(macroTotal)
        xs = self % total

      case(macroCapture)
        xs = self % capture

      case(macroEscatter)
        xs = self % elasticScatter

      !case(macroAbsorbtion)
      !  xs = self % fission + self % capture

      case default
        xs = ZERO

    end select

  end function get

  !!
  !! Use a real r in <0;1> to sample reaction from Macroscopic XSs
  !!
  !! This function might be common thus is type-bound procedure for conveniance
  !!
  !! Args:
  !!   r [in] -> Real number in <1;0>
  !!
  !! Result:
  !!   One of the Macroscopic MT numbers
  !!     elasticScatter   = macroEscatter
  !!     inelasticScatter = macroIEscatter
  !!     capture          = macroCapture
  !!
  !! Errors::
  !!   If r < 0 then returns macroEscatter
  !!
  elemental function invert_macroXSs(self, r) result(MT)
    class(IMCMacroXSs), intent(in) :: self
    real(defReal), intent(in)          :: r
    integer(shortInt)                  :: MT
    real(defReal)                      :: xs
    integer(shortInt)                  :: C

    ! Elastic Scattering
    C = 1
    xs = self % total * r - self % elasticScatter
    if (xs > ZERO) C = C + 1

    ! Inelastic Scattering
    xs = xs - self % inelasticScatter
    if(xs > ZERO) C = C + 1

    ! Capture
    xs = xs - self % capture
    if(xs > ZERO) C = C + 1

    ! Choose MT number
    select case(C)
      case(1)
        MT = macroEScatter

      case(2)
        MT = macroIEscatter

      case(3)
        MT = macroCapture

      case default  ! Should never happen -> Avoid compiler error and return nonsense number
        MT = huge(C)

    end select

  end function invert_macroXSs


  !!
  !! Use a real r in <0;1> to sample reaction from Microscopic XSs
  !!
  !! This function involves a bit of code so is written for conviniance
  !!
  !! Args:
  !!   r [in] -> Real number in <0;1>
  !!
  !! Result:
  !!   MT number of the reaction:
  !!     elastic scatter   = N_N_elastic
  !!     inelastic scatter = N_N_inelastic
  !!     capture           = N_diasp
  !!
  !! Errors:
  !!   If r < 0 then returns N_N_elastic
  !!
  elemental function invert_microXSs(self, r) result(MT)
    class(IMCMicroXSs), intent(in) :: self
    real(defReal), intent(in)          :: r
    integer(shortInt)                  :: MT
    real(defReal)                      :: xs
    integer(shortInt)                  :: C

    ! Elastic Scattering
    C = 1
    xs = self % total * r - self % elasticScatter
    if (xs > ZERO) C = C + 1

    ! Inelastic Scattering
    xs = xs - self % inelasticScatter
    if(xs > ZERO) C = C + 1

    ! Capture
    xs = xs - self % capture
    if(xs > ZERO) C = C + 1

    ! Choose MT number
    select case(C)
      case(1)
        MT = N_N_elastic

      case(2)
        MT = N_N_inelastic

      case(3)
        MT = N_disap

      case default  ! Should never happen -> Avoid compiler error and return nonsense number
        MT = huge(C)
    end select

  end function invert_microXSs


end module IMCXsPackages_class
