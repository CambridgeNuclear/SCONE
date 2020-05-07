module fluxResponse_class

  use numPrecision
  use dictionary_class,    only : dictionary
  use particle_class,      only : particle
  use tallyResponse_inter, only : tallyResponse

  ! Nuclear Data interface
  use nuclearDatabase_inter, only : nuclearDatabase

  implicit none
  private

  !!
  !! tallyResponse to score flux contribution
  !!
  !! Always returns ONE
  !!
  !! Interface:
  !!   tallyResponse Interface
  !!
  type, public,extends(tallyResponse) :: fluxResponse
    private
  contains
    procedure :: init
    procedure :: get
    procedure :: kill
  end type fluxResponse

contains

  !!
  !! Initialise Response from dictionary
  !!
  !! See tallyResponse_inter for details
  !!
  subroutine init(self, dict)
    class(fluxResponse), intent(inout) :: self
    class(dictionary), intent(in)      :: dict

    ! Do nothing

  end subroutine init

  !!
  !! Get 1.0 (Response to score flux)
  !!
  !! See tallyResponse_inter for details
  !!
  function get(self, p, xsData) result(val)
    class(fluxResponse), intent(in)       :: self
    class(particle), intent(in)           :: p
    class(nuclearDatabase), intent(inout) :: xsData
    real(defReal)                         :: val

    val = ONE

  end function get

  !!
  !! Return to uninitialised State
  !!
  elemental subroutine kill(self)
    class(fluxResponse), intent(inout) :: self

    ! Do nothing for nothing can be done

  end subroutine kill

end module fluxResponse_class
