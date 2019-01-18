module fluxResponse_class

  use numPrecision
  use dictionary_class,    only : dictionary
  use particle_class,      only : particle
  use tallyResponse_inter, only : tallyResponse

  implicit none
  private

  !!
  !! tallyResponse to score flux contribution
  !!  Always returns ONE
  !!  Should be used only for testing really...
  !!
  type, public,extends(tallyResponse) :: fluxResponse
    private
  contains
    procedure  :: init
    procedure  :: get
  end type fluxResponse

contains

  !!
  !! Initialise Response from dictionary
  !!
  subroutine init(self, dict)
    class(fluxResponse), intent(inout) :: self
    class(dictionary), intent(in)      :: dict

    ! Do nothing

  end subroutine init

  !!
  !! Get 1.0 (Response to score flux)
  !!
  elemental function get(self, p) result(val)
    class(fluxResponse), intent(in) :: self
    class(particle), intent(in)     :: p
    real(defReal)                   :: val

    val = ONE

  end function get

end module fluxResponse_class
