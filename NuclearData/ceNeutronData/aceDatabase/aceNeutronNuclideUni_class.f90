module aceNeutronNuclideUni_class

  use numPrecision

  ! Nuclear Data Interfaces
  use aceNeutronNuclide_class,      only : aceNeutronNuclide
  use neutronXSPackages_class,      only : neutronMicroXSs


  implicit none
  private


  ! Grid location parameters
  integer(shortInt), parameter :: TOTAL_XS      = 1
  integer(shortInt), parameter :: ESCATTER_XS   = 2
  integer(shortInt), parameter :: IESCATTER_XS  = 3
  integer(shortInt), parameter :: CAPTURE_XS    = 4
  integer(shortInt), parameter :: FISSION_XS    = 5
  integer(shortInt), parameter :: NU_FISSION    = 6


  !!
  !! Extension of aceNeutronNuclide, used to read and write material unionised energy grids.
  !! For more information see aceNeutronNuclide_class
  !!
  !! Called only by aceNeutronDatabaseUni
  !!
  !! Public Members:
  !!   unionData       -> Array of XSs that are required in ceNeutronMicroXSs, that is
  !!     (total, capture, escatter, iescatter, fission, nuFission)
  !!   matIdx          -> Material the unionised grid referes to. Used by aceNeutronDatabaseUni
  !!
  !! Procedures:
  !!   totalXS  -> return totalXS given index and interpolation factor
  !!   microXSs -> return interpolated ceNeutronMicroXSs package given index and inter. factor
  !!
  type, public, extends(aceNeutronNuclide) :: aceNeutronNuclideUni
    real(defReal), dimension(:,:), allocatable  :: unionData
    integer(shortInt)                           :: matIdx  = 0
  contains
    ! Override superclass procedures
    procedure :: totalXS
    procedure :: microXSs
  end type aceNeutronNuclideUni

contains

  !!
  !! Return value of the total XS given interpolation factor and index
  !!
  !! Does not prefeorm any check for valid input!
  !!
  !! For more information see aceNeutronNuclide_class
  !!
  elemental function totalXS(self, idx, f) result(xs)
    class(aceNeutronNuclideUni), intent(in) :: self
    integer(shortInt), intent(in)           :: idx
    real(defReal), intent(in)               :: f
    real(defReal)                           :: xs

    xs = self % unionData(TOTAL_XS, idx+1) * f + (ONE-f) * self % unionData(TOTAL_XS, idx)

  end function totalXS

  !!
  !! Return interpolated neutronMicroXSs package for the given interpolation factor and index
  !!
  !! Does not prefeorm any check for valid input!
  !!
  !! For more information see aceNeutronNuclide_class
  !!
  elemental subroutine microXSs(self, xss, idx, f)
    class(aceNeutronNuclideUni), intent(in) :: self
    type(neutronMicroXSs), intent(out)      :: xss
    integer(shortInt), intent(in)           :: idx
    real(defReal), intent(in)               :: f

    associate (data => self % unionData(:,idx:idx+1))

      xss % total            = data(TOTAL_XS, 2)  * f + (ONE-f) * data(TOTAL_XS, 1)
      xss % elasticScatter   = data(ESCATTER_XS, 2)  * f + (ONE-f) * data(ESCATTER_XS, 1)
      xss % inelasticScatter = data(IESCATTER_XS, 2) * f + (ONE-f) * data(IESCATTER_XS, 1)
      xss % capture          = data(CAPTURE_XS, 2)   * f + (ONE-f) * data(CAPTURE_XS, 1)

      if(self % isFissile()) then
        xss % fission   = data(FISSION_XS, 2) * f + (ONE-f) * data(FISSION_XS, 1)
        xss % nuFission = data(NU_FISSION, 2) * f + (ONE-f) * data(NU_FISSION, 1)
      else
        xss % fission   = ZERO
        xss % nuFission = ZERO
      end if
    end associate

  end subroutine microXSs

end module aceNeutronNuclideUni_class
