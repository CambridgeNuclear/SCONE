module testResponse_class

  use numPrecision
  use dictionary_class,    only : dictionary
  use particle_class,      only : particle
  use tallyResponse_inter, only : tallyResponse

  implicit none
  private

  !!
  !! tallyResponse to be used for testing
  !!  Always returns a constant
  !!
  type, public,extends(tallyResponse) :: testResponse
    private
    real(defReal) :: value = ONE
  contains
    procedure  :: init
    procedure  :: get
  end type testResponse

contains

  !!
  !! Initialise Response from dictionary
  !!
  subroutine init(self, dict)
    class(testResponse), intent(inout) :: self
    class(dictionary), intent(in)      :: dict

    call dict % get(self % value,'value')

  end subroutine init

  !!
  !! Get 1.0 (Response to score flux)
  !!
  function get(self, p) result(val)
    class(testResponse), intent(in) :: self
    class(particle), intent(in)     :: p
    real(defReal)                   :: val

    val = self % value

  end function get

end module testResponse_class
