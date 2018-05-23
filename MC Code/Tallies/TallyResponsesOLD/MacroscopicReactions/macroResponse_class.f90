module macroResponse_class

  use numPrecision
  use endfConstants
  use genericProcedures,          only : fatalError
  use dictionary_class,           only : dictionary
  use tallyResponse_inter,    only : tallyResponse
  use particle_class,             only : particle, phaseCoord
  use nuclearData_inter,          only : nuclearData
  use transportNuclearData_inter, only : transportNuclearData

  use xsMacroSet_class,           only : xsMacroSet_ptr

  implicit none
  private

  type, public,extends(tallyResponse) :: macroResponse
    private
  contains
    procedure :: init
    procedure :: getScore
  end type macroResponse

contains

  !!
  !! Response Constructor
  !!
  function macroResponse_cont(dict) result(response)
    class(dictionary), intent(in)     :: dict
    class(tallyResponse),allocatable  :: response
    type(macroResponse),allocatable   :: locResponse

    allocate( locResponse)
    call locResponse % init(dict)

    call move_alloc(locResponse, response)

  end function macroResponse_cont

  !!
  !! Initialise from dictionary
  !!
  subroutine init(self,dict)
    class(macroResponse),intent(inout) :: self
    class(dictionary),intent(in)       :: dict

    ! Get response code
    call dict % get(self % responseCode,'MT')

  end subroutine init

  !!
  !! Response multiplier of flux (f*flux)
  !! Is equal to macroscopic XS in collision material at collison energy
  !!
  function getScore(self,pre,post,MT,muL) result (f)
    class(macroResponse), intent(inout) :: self
    class(phaseCoord), intent(in)       :: pre
    class(particle), intent(in)         :: post
    integer(shortInt), intent(in)       :: MT
    real(defReal), intent(in)           :: muL
    real(defReal)                       :: f
    type(xsMacroSet_ptr)                :: macroXss
    character(100), parameter           :: Here =' getCollisionScore (macroResponse_class.f90)'

    ! Obtain pointer to data from particle
    ! Check if it dynamic type is supported
    ! If it is obtain macroscopic XSs
    ! It it isn't throw error
    associate (xsData => post % xsData)
      select type(xsData)
        class is (transportNuclearData)
          call xsData % getMatMacroXS(macroXss, post, post % matIdx)

        class default
          call fatalError(Here,'Dynamic type of XS data attached to particle is not transportNuclearData')

      end select
    end associate

    ! Select approperiate Macroscopic XS based on response code
    select case(self % responseCode)
      case(macroTotal)
        f = macroXss % totalXS()

      case(macroCapture)
        f = macroXss % captureXS()

      case(macroAllScatter)
        f = macroXss % scatterXS()

      case(macroAbsorbtion)
        f = macroXss % captureXS() + macroXss % fissionXS()

      case(macroFission)
        f = macroXss % fissionXS()

      case(macroNuFission)
        f = macroXss % nuFissionXS()

    end select


  end function getScore
    
end module macroResponse_class
