module macroResponse_class

  use numPrecision
  use endfConstants
  use universalVariables,         only : VOID_MAT
  use genericProcedures,          only : fatalError, numToChar
  use dictionary_class,           only : dictionary
  use particle_class,             only : particle, P_NEUTRON
  use tallyResponse_inter,        only : tallyResponse

  ! Nuclear Data interfaces
  use nuclearDatabase_inter,      only : nuclearDatabase
  use neutronMaterial_inter,      only : neutronMaterial, neutronMaterial_CptrCast
  use neutronXsPackages_class,    only : neutronMacroXSs

  implicit none
  private


  !!
  !! tallyResponse for scoring a single macroscopicXSs
  !!  Currently supports neutrons only
  !!
  !! Private Members:
  !!   MT -> MT number of the macroscopic reaction for weighting
  !!
  !! Interface:
  !!   tallyResponse interface
  !!   build -> Initialise directly from MT number
  !!
  !! Sample dictionary input
  !!  name {
  !!     type macroResponse;
  !!     MT   <int>;
  !!  }
  !!
  type, public,extends(tallyResponse) :: macroResponse
    private
    !! Response MT number
    integer(shortInt) :: MT = 0
  contains
    ! Superclass Procedures
    procedure  :: init
    procedure  :: get
    procedure  :: kill

    ! Local Procedures
    procedure  :: build

  end type macroResponse

contains

  !!
  !! Initialise Response from dictionary
  !!
  !! See tallyResponse_inter for details
  !!
  !! Errors:
  !!   fatalError if MT is invalid
  !!
  subroutine init(self, dict)
    class(macroResponse), intent(inout) :: self
    class(dictionary), intent(in)       :: dict
    integer(shortInt)                   :: MT
    character(100), parameter :: Here = 'init ( macroResponse_class.f90)'

    ! Load MT number
    call dict % get(MT, 'MT')

    ! Build response
    call self % build(MT)

  end subroutine init

  !!
  !! Build macroResponse from MT number
  !!
  !! Args:
  !!   MT [in] -> MT number for weighting
  !!
  !! Errors:
  !!   fatalError if MT is invalid
  !!
  subroutine build(self, MT)
    class(macroResponse), intent(inout) :: self
    integer(shortInt), intent(in)       :: MT
    character(100), parameter :: Here = 'build ( macroResponse_class.f90)'

    ! Check that MT number is valid
    select case(MT)
      case(macroTotal, macroCapture, macroFission, macroNuFission, macroAbsorbtion)
        ! Do nothing. MT is Valid

      case(macroEscatter)
        call fatalError(Here,'Macroscopic Elastic scattering is not implemented yet')

      case default
        call fatalError(Here,'Unrecognised MT number: '// numToChar(self % MT))
    end select

    ! Load MT
    self % MT = MT

  end subroutine build

  !!
  !! Return response value
  !!
  !! See tallyResponse_inter for details
  !!
  !! Errors:
  !!   Return ZERO if particle is not a Neutron
  !!
  function get(self, p, xsData) result(val)
    class(macroResponse), intent(in)      :: self
    class(particle), intent(in)           :: p
    class(nuclearDatabase), intent(inout) :: xsData
    real(defReal)                         :: val
    type(neutronMacroXSs)                 :: xss
    class(neutronMaterial), pointer       :: mat

    val = ZERO

    ! Return zero if particle is not neutron or if the particle is in void
    if (p % type /= P_NEUTRON) return
    if (p % matIdx() == VOID_MAT) return

    ! Get pointer to active material data
    mat => neutronMaterial_CptrCast(xsData % getMaterial(p % matIdx()))

    ! Return if material is not a neutronMaterial
    if (.not.associated(mat)) return

    call mat % getMacroXSs(xss, p)
    val = xss % get(self % MT)

  end function get

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(macroResponse), intent(inout) :: self

    self % MT = 0

  end subroutine kill

end module macroResponse_class
