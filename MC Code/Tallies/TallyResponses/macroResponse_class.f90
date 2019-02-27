module macroResponse_class

  use numPrecision
  use endfConstants
  use genericProcedures,          only : fatalError, numToChar
  use dictionary_class,           only : dictionary
  use particle_class,             only : particle
  use tallyResponse_inter,        only : tallyResponse

  ! Nuclear Data interface
  use transportNuclearData_inter, only : transportNuclearData
  use xsMacroSet_class,           only : xsMacroSet_ptr


  implicit none
  private


  !!
  !! tallyResponse for scoring a single macroscopicXSs
  !!  Currently supports neutrons only
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
    procedure  :: init
    procedure  :: get
  end type macroResponse

contains

  !!
  !! Initialise Response from dictionary
  !!
  subroutine init(self, dict)
    class(macroResponse), intent(inout) :: self
    class(dictionary), intent(in)       :: dict
    character(100), parameter :: Here = 'init ( macroResponse_class.f90)'

    ! Load MT number
    call dict % get(self % MT, 'MT')

    ! Check that MT number is valid
    select case(self % MT)
      case(macroTotal, macroCapture, macroFission, macroNuFission, macroAbsorbtion)
        ! Do nothing. MT is Valid

      case(macroEscatter)
        call fatalError(Here,'Macroscopic Elastic scattering is not implemented yet')

      case default
        call fatalError(Here,'Unrecognised MT number: '// numToChar(self % MT))
    end select


  end subroutine init

  !!
  !! Return response value
  !!  Returns 0.0 if the xs type is invalid
  !!
  function get(self, p) result(val)
    class(macroResponse), intent(in) :: self
    class(particle), intent(in)      :: p
    real(defReal)                    :: val
    type(xsMacroSet_ptr)             :: macroXSs

    select type(xs => p % xsData)
      class is(transportNuclearData)
        call xs % getMatMacroXS(macroXSs,p , p % matIdx())
        val = macroXSs % xsOf(self % MT)

      class default
        val = ZERO

    end select
  end function get
    
end module macroResponse_class
