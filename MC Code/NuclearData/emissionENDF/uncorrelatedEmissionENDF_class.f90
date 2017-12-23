module uncorrelatedEmissionENDF_class

  use numPrecision
  use emissionENDF_class, only : emissionENDF
  use RNG_class, only :RNG

  implicit none
  private

  type, public, extends(emissionENDF) :: uncorrelatedEmissionENDF
      private
      ! Pointer to EnergyENDFLaw
      ! Pointer to AngleENDFLaw
      ! Pointer to neutronReleseENDF
    contains
      procedure :: getAngleEnergy
      procedure :: getNumber
      ! Build procedures
      generic   :: attachENDF   => attachENDF_Angle , &
                                   attachENDF_Energy, &
                                   attachENDF_Relese
      ! Private procedures
      procedure,private :: attachENDF_Angle
      procedure,private :: attachENDF_Energy
      procedure,private :: attachENDF_Relese


  end type uncorrelatedEmissionENDF

contains

  subroutine getAngleEnergy(self,angle,energy,rand )
    !! Subroutine, which returns a sample of angle and energy obtained from law attached to the
    !! object.
    class(uncorrelatedEmissionENDF), intent(in)  :: self
    real(defReal), intent(inout)                 :: angle
    real(defReal), intent(inout)                 :: energy
    class(RNG), intent(inout)                    :: rand
  end subroutine

  subroutine getNumber(self,number)
    !! Subroutine, which returns a number of emitted secondary neutrons according to the attached
    !! neutronReleseENDF object.
    class(uncorrelatedEmissionENDF), intent(in) :: self
    real(defReal), intent(inout)                :: number
  end subroutine

  subroutine attachENDF_Angle(self,A)
    !! Subroutine, which attaches pointer to angular distribution of emmited neutrons
    class(uncorrelatedEmissionENDF), intent(inout) :: self
    integer :: A
  end subroutine attachENDF_Angle

   subroutine attachENDF_Energy(self,E)
    !! Subroutine, which attaches pointer to energy distribution of emmited neutrons
    class(uncorrelatedEmissionENDF), intent(inout) :: self
    real :: E
  end subroutine attachENDF_Energy

  subroutine attachENDF_Relese(self,B)
    !! Subroutine which attaches pointer to distribution of secondary neutrons
    class(uncorrelatedEmissionENDF), intent(inout) :: self
    logical :: B
  end subroutine attachENDF_Relese

end module uncorrelatedEmissionENDF_class
