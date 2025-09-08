module invSpeedResponse_class

  use numPrecision
  use universalVariables,  only : neutronMass, lightSpeed
  use dictionary_class,    only : dictionary
  use particle_class,      only : particle, P_NEUTRON, P_PHOTON
  use tallyResponse_inter, only : tallyResponse

  ! Nuclear Data interface
  use nuclearDatabase_inter, only : nuclearDatabase

  implicit none
  private

  !!
  !! tallyResponse which returns the inverse of the particle speed in [cm/s]
  !!
  !! Can be also used to estimate particle density
  !!
  !! NOTE:
  !!  The speeds are computed from non-relativistic formula for massive particles.
  !!  A small error might appear in MeV range (e.g. for fusion applications)
  !!
  !! Interface:
  !!   tallyResponse Interface
  !!
  type, public,extends(tallyResponse) :: invSpeedResponse
    private
  contains
    procedure :: init
    procedure :: get
    procedure :: kill

  end type invSpeedResponse

contains

  !!
  !! Initialise Response from dictionary
  !!
  !! See tallyResponse_inter for details
  !!
  subroutine init(self, dict)
    class(invSpeedResponse), intent(inout) :: self
    class(dictionary), intent(in)         :: dict

    ! Do nothing

  end subroutine init

  !!
  !! Returns the inverse of the particle speed
  !!
  !! See tallyResponse_inter for details
  !!
  function get(self, p, xsData) result(val)
    class(invSpeedResponse), intent(in)    :: self
    class(particle), intent(in)           :: p
    class(nuclearDatabase), intent(inout) :: xsData
    real(defReal)                         :: val

    ! Gets the particle speed from the particle
    val = ONE / p % getSpeed()

  end function get

  !!
  !! Return to uninitialised State
  !!
  elemental subroutine kill(self)
    class(invSpeedResponse), intent(inout) :: self

    ! Do nothing for nothing can be done

  end subroutine kill

end module invSpeedResponse_class
