module weightResponse_class

  use numPrecision
  use endfConstants
  use genericProcedures,          only : fatalError, numToChar
  use dictionary_class,           only : dictionary
  use particle_class,             only : particle, P_NEUTRON
  use tallyResponse_inter,        only : tallyResponse

  ! Nuclear Data interfaces
  use nuclearDatabase_inter,      only : nuclearDatabase
  use neutronMaterial_inter,      only : neutronMaterial, neutronMaterial_CptrCast

  implicit none
  private


  !!
  !! tallyResponse for scoring particle weights
  !!  Currently supports neutrons only
  !!
  !! Interface:
  !!   tallyResponse interface
  !!
  !! Sample dictionary input
  !!  name {
  !!     type weightResponse; moment 1;
  !!  }
  !!
  type, public,extends(tallyResponse) :: weightResponse
    private
    integer(shortInt)    :: moment
  contains
    ! Superclass Procedures
    procedure  :: init
    procedure  :: get
    procedure  :: kill

  end type weightResponse

contains

  !!
  !! Initialise Response from dictionary
  !!
  !! See tallyResponse_inter for details
  !!
  subroutine init(self, dict)
    class(weightResponse), intent(inout) :: self
    class(dictionary), intent(in)        :: dict
    character(100), parameter :: Here ='init (weightResponse_class.f90)'

    ! Get response moment to be calculated
    call dict % getOrDefault(self % moment, 'moment', 1)

    if (self % moment < 0) call fatalError(Here, 'Moment must be bigger or equal zero.')

  end subroutine init

  !!
  !! Return response value
  !!
  !! See tallyResponse_inter for details
  !!
  !! Errors:
  !!   Return ZERO if particle is not a Neutron
  !!
  function get(self, p, xsData) result(val)
    class(weightResponse), intent(in)      :: self
    class(particle), intent(in)           :: p
    class(nuclearDatabase), intent(inout) :: xsData
    real(defReal)                         :: val
    class(neutronMaterial), pointer       :: mat

    val = ZERO

    ! Return 0.0 if particle is not neutron
    if(p % type /= P_NEUTRON) return

    ! Get pointer to active material data
    mat => neutronMaterial_CptrCast(xsData % getMaterial(p % matIdx()))

    ! Return if material is not a neutronMaterial
    if(.not.associated(mat)) return

    if (self % moment == 0) then
      val = xsData % getTotalMatXS(p, p % matIdx()) / (p % w)
    else
      val = xsData % getTotalMatXS(p, p % matIdx()) * ((p % w) ** (self % moment - 1))
    end if

  end function get

  !!
  !! Return to uninitialised State
  !!
  elemental subroutine kill(self)
    class(weightResponse), intent(inout) :: self

    ! Do nothing for nothing can be done

  end subroutine kill

end module weightResponse_class
