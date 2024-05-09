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
  !! Returns the inverse of the particle velocity in [cm/s]
  !!
  !! The velocity is calculated from the particle energy for neutrons, and it is
  !! the speed of light for photons
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
  !! Calculates the particle velocity for neutron and photons
  !! NOTE: neutronMass: [MeV]
  !!       lightSpeed:  [cm/s]
  !!
  !! The function returns the inverse of the velocity (response to score particle density)
  !! if the particle type is neutron or photon, and zero otherwise
  !!
  !! See tallyResponse_inter for details
  !!
  function get(self, p, xsData) result(val)
    class(densityResponse), intent(in)    :: self
    class(particle), intent(in)           :: p
    class(nuclearDatabase), intent(inout) :: xsData
    real(defReal)                         :: val
    real(defReal)                         :: velocity

    ! Initialise response
    val = ZERO

    ! Calculates the velocity for the relevant particle [cm/s]
    if (p % type == P_NEUTRON) then
      velocity = sqrt(TWO * p % E / neutronMass) * lightSpeed

    elseif (p % type == P_PHOTON) then
      velocity = lightSpeed

    else
      return

    end if

    ! Returns the inverse of the velocity
    val = ONE / velocity

  end function get

  !!
  !! Return to uninitialised State
  !!
  elemental subroutine kill(self)
    class(densityResponse), intent(inout) :: self

    ! Do nothing for nothing can be done

  end subroutine kill

end module densityResponse_class
