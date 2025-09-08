module microResponse_class

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

  ! Material interface
  use materialMenu_mod,           only : materialItem, matName, nMat, getMatPtr

  implicit none
  private


  !!
  !! tallyResponse for scoring a single microscopicXSs
  !!   Currently supports neutrons only
  !!
  !! Private Members:
  !!   MT     -> MT number of the microscopic reaction for weighting
  !!   matIdx -> index of the material that contains (only) the nuclide wanted
  !!   dens   -> atomic density of the nuclide
  !!
  !! Interface:
  !!   tallyResponse interface
  !!   build -> Initialise directly from MT number
  !!
  !! Sample dictionary input:
  !!  name {
  !!     type microResponse;
  !!     MT   <int>;
  !!     material <matName>;
  !!  }
  !!
  !! Note:
  !!   The material <matName> must include only one nuclide. The final estimate
  !!   is independent of the nuclide atomic density, which can be any value but zero.
  !!
  type, public,extends(tallyResponse) :: microResponse
    private
    !! Response MT number
    integer(shortInt)        :: MT
    integer(shortInt),public :: matIdx
    real(defReal)            :: dens
  contains
    ! Superclass Procedures
    procedure  :: init
    procedure  :: get
    procedure  :: kill

    ! Local Procedures
    procedure  :: build

  end type microResponse

contains

  !!
  !! Initialise Response from dictionary
  !!
  !! See tallyResponse_inter for details
  !!
  !! Errors:
  !!   fatalError if the material contains more than one nuclide
  !!   fatalError if the nuclide has density 0.0
  !!
  subroutine init(self, dict)
    class(microResponse), intent(inout) :: self
    class(dictionary), intent(in)       :: dict
    integer(shortInt)                   :: MT, i
    character(15)                       :: mName
    type(materialItem), pointer         :: mat
    character(100), parameter :: Here = 'init ( microResponse_class.f90)'

    ! Load MT number and material name
    call dict % get(MT, 'MT')
    call dict % get(mName, 'material')

    ! Find corresponding material index
    do i = 1,nMat()
      if (mName == matName(i)) self % matIdx = i
    end do

    ! Get pointer to the material
    mat => getMatPtr(self % matIdx)

    if (size(mat % dens) > 1) call fatalError(Here, 'Material '//trim(mName)//' &
                                              & has more than one nuclide' )
    self % dens = mat % dens(1)

    if (self % dens == ZERO) call fatalError(Here, 'Density of material &
                                             & '//trim(mName)//' cannot be 0' )
    ! Build response
    call self % build(MT)

  end subroutine init

  !!
  !! Build microResponse from MT number
  !!
  !! Args:
  !!   MT [in] -> MT number for weighting
  !!
  !! Errors:
  !!   fatalError if MT is invalid
  !!
  subroutine build(self, MT)
    class(microResponse), intent(inout) :: self
    integer(shortInt), intent(in)       :: MT
    character(100), parameter :: Here = 'build ( microResponse_class.f90)'

    ! Check that MT number is valid and load MT
    select case(MT)
      case(N_TOTAL)
        self % MT = macroTotal
      case(N_N_ELASTIC)
        self % MT = macroEScatter
      case(N_GAMMA)
        self % MT = macroCapture
      case(N_FISSION)
        self % MT = macroFission
      case(N_ABSORPTION)
        self % MT = macroAbsorbtion
      case default
        call fatalError(Here,'Unrecognised MT number: '// numToChar(MT))
    end select

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
    class(microResponse), intent(in)      :: self
    class(particle), intent(in)           :: p
    class(nuclearDatabase), intent(inout) :: xsData
    real(defReal)                         :: val
    class(neutronMaterial), pointer       :: mat
    type(neutronMacroXSs)                 :: xss
    character(100), parameter :: Here = 'get ( microResponse_class.f90)'

    val = ZERO

    ! Return zero if particle is not neutron or if the particle is in void
    if (p % type /= P_NEUTRON) return
    if (p % matIdx() == VOID_MAT) return

    ! Get pointer to active material data
    mat => neutronMaterial_CptrCast(xsData % getMaterial(self % matIdx))

    ! Return if material is not a neutronMaterial
    if (.not.associated(mat)) return

    ! Get the macroscopic cross section for the material
    call mat % getMacroXSs(xss, p)

    ! Normalise the macroscopic cross section with the atomic density
    val = xss % get(self % MT) / self % dens

  end function get

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(microResponse), intent(inout) :: self

    self % MT = 0

  end subroutine kill

end module microResponse_class
