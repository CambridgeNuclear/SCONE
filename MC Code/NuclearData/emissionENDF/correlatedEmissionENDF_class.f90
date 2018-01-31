module correlatedEmissionENDF_class

  use numPrecision
  use emissionENDF_class, only : emissionENDF
  use RNG_class, only :RNG

  implicit none
  private

  type, public, extends(emissionENDF) :: correlatedEmissionENDF
      private
      ! Pointer to correletedENDFLaw
      ! Pointer to neutronReleseENDF
    contains
      procedure :: sampleAngleEnergy
      procedure :: releaseAt
      ! Build procedures
      generic   :: attachENDF => attachENDF_Correlated, &
                                 attachENDF_Relese
      ! Private procedures
      procedure,private :: attachENDF_Correlated
      procedure,private :: attachENDF_Relese

  end type correlatedEmissionENDF

contains
  subroutine sampleAngleEnergy(self,angle,E_out,E_in,rand )
    !! Subroutine, which returns a sample of angle and energy obtained from law attached to the
    !! object.
    class(correlatedEmissionENDF), intent(in)    :: self
    real(defReal), intent(out)                   :: angle
    real(defReal), intent(out)                   :: E_out
    real(defReal), intent(in)                    :: E_in
    class(RNG), intent(inout)                    :: rand

  end subroutine

  function releaseAt(self,E_in) result(number)
    !! Subroutine, which returns a number of emitted secondary neutrons according to the attached
    !! neutronReleseENDF object.
    class(correlatedEmissionENDF), intent(in)   :: self
    real(defReal), intent(in)                   :: E_in
    real(defReal)                               :: number
  end function releaseAt

  subroutine attachENDF_Correlated(self,A)
    !! Subroutine, which attaches pointer to angular distribution of emmited neutrons
    class(correlatedEmissionENDF), intent(inout) :: self
    integer :: A
  end subroutine attachENDF_Correlated

  subroutine attachENDF_Relese(self,B)
    !! Subroutine which attaches pointer to distribution of secondary neutrons
    class(correlatedEmissionENDF), intent(inout) :: self
    logical :: B
  end subroutine attachENDF_Relese
    
end module correlatedEmissionENDF_class
