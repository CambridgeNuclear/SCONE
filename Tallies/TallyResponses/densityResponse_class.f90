module densityResponse_class

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
  !! tallyResponse to score particle density contribution
  !!
  !! Returns the inverse of the particle speed in [cm/s]
  !!
  !! NOTE:
  !!  The speeds are computed from non-relativistic formula for massive particles.
  !!  The small error might appear in MeV range (e.g. for fusion applications)
  !!
  !! Interface:
  !!   tallyResponse Interface
  !!
  type, public,extends(tallyResponse) :: densityResponse
    private
  contains
    procedure :: init
    procedure :: get
    procedure :: kill

  end type densityResponse

contains

  !!
  !! Initialise Response from dictionary
  !!
  !! See tallyResponse_inter for details
  !!
  subroutine init(self, dict)
    class(densityResponse), intent(inout) :: self
    class(dictionary), intent(in)         :: dict

    ! Do nothing

  end subroutine init

  !!
  !! Returns the inverse of the particle speed (response to score particle density)
  !!
  !! See tallyResponse_inter for details
  !!
  function get(self, p, xsData) result(val)
    class(densityResponse), intent(in)    :: self
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
    class(densityResponse), intent(inout) :: self

    ! Do nothing for nothing can be done

  end subroutine kill

end module densityResponse_class
