module uncorrelatedEmissionENDF_class

  use numPrecision
  use emissionENDF_class,   only : emissionENDF
  use RNG_class,            only : RNG
  use energyLawENDF_class,  only : energyLawENDF
  use angleLawENDF_class,   only : angleLawENDF
  use releaseLawENDF_class, only : releaseLawENDF

  implicit none
  private

  interface uncorrelatedEmissionENDF
    module procedure new_uncorrelatedEmissionENDF
  end interface


  type, public, extends(emissionENDF) :: uncorrelatedEmissionENDF
      private

      class(angleLawENDF), pointer   :: anglePdf   => null()
      class(energyLawENDF), pointer  :: energyPdf  => null()
      class(releaseLawENDF), pointer :: releasePdf => null() !! Not exactly a PDF (its a delta function PDF...)

    contains
      procedure :: getAngleEnergy
      procedure :: getNumber
      ! Build procedures
      generic   :: attachENDF   => attachENDF_Angle , &
                                   attachENDF_Energy, &
                                   attachENDF_Release
      ! Private procedures
      procedure,private :: attachENDF_Angle
      procedure,private :: attachENDF_Energy
      procedure,private :: attachENDF_Release


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


  subroutine attachENDF_Angle(self,anglePdf)
    !! Subroutine, which attaches pointer to angular distribution of emmited neutrons
    class(uncorrelatedEmissionENDF), intent(inout) :: self
    class(angleLawENDF),target, intent(in)         :: anglePdf

    if(associated(self % anglePdf)) deallocate(self % anglePdf)

    self % anglePdf => anglePdf

  end subroutine attachENDF_Angle


  subroutine attachENDF_Energy(self,energyPdf)
    !! Subroutine, which attaches pointer to energy distribution of emmited neutrons
    class(uncorrelatedEmissionENDF), intent(inout)  :: self
    class(energyLawENDF),target, intent(in)         :: energyPdf

    if(associated(self % energyPdf)) deallocate(self % energyPdf)

    self % energyPdf => energyPdf

  end subroutine attachENDF_Energy


  subroutine attachENDF_Release(self,releasePdf)
    !! Subroutine which attaches pointer to distribution of secondary neutrons
    class(uncorrelatedEmissionENDF), intent(inout) :: self
    class(releaseLawENDF),target, intent(in)       :: releasePdf

    if(associated(self % releasePdf)) deallocate(self % releasePdf)

    self % releasePdf => releasePdf

  end subroutine attachENDF_Release


  function new_uncorrelatedEmissionENDF(angleLaw,energyLaw,releaseLaw) result (new)
    class(angleLawENDF),pointer, intent(in)       :: angleLaw
    class(energyLawENDF),pointer,intent(in)       :: energyLaw
    class(releaseLawENDF),pointer,intent(in)      :: releaseLaw
    type(uncorrelatedEmissionENDF),pointer        :: new

    allocate(new)

    call new % attachENDF(angleLaw)
    call new % attachENDF(energyLaw)
    call new % attachENDF(releaseLaw)

  end function new_uncorrelatedEmissionENDF


end module uncorrelatedEmissionENDF_class
