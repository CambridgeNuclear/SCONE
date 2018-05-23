module fluxResponse_class

  use numPrecision
  use endfConstants
  use dictionary_class,    only : dictionary
  use tallyResponse_inter, only : tallyResponse
  use particle_class,      only : particle, phaseCoord

  implicit none
  private

  type, public,extends(tallyResponse) :: fluxResponse
    private
  contains
    procedure :: init
    procedure :: getScore
  end type fluxResponse

contains

  !!
  !! Response Constructor
  !!
  function fluxResponse_cont(dict) result(response)
    class(dictionary), intent(in) :: dict
    class(tallyResponse),allocatable :: response
    type(fluxResponse),allocatable   :: locResponse

    allocate( locResponse)
    call locResponse % init(dict)

    call move_alloc(locResponse, response)

  end function fluxResponse_cont

  !!
  !! Initialisation procedure
  !! In this case it virtually nothing
  !!
  subroutine init(self,dict)
    class(fluxResponse), intent(inout) :: self
    class(dictionary), intent(in)      :: dict

    self % responseCode = macroTotal

  end subroutine init

  !!
  !! Response multiplier of flux (f*flux)
  !! To get flux it needs to be ONE
  !! Very underwhelming function
  !!
  function getScore(self,pre,post,MT,muL) result (f)
    class(fluxResponse), intent(inout)  :: self
    class(phaseCoord), intent(in)       :: pre
    class(particle), intent(in)         :: post
    integer(shortInt), intent(in)       :: MT
    real(defReal), intent(in)           :: muL
    real(defReal)                       :: f

    f = ONE

  end function getScore
    
end module fluxResponse_class
