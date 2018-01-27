module emissionENDF_class

  use numPrecision
  use RNG_class

  implicit none
  private

  type,abstract, public :: emissionENDF
    private
    integer(shortInt)  :: MT
    logical(defBool)   :: cmFrame = .true.
  contains
    procedure(getAngleEnergy),deferred :: getAngleEnergy
    procedure(getNumber),deferred      :: getNumber
    procedure                          :: setMT
    procedure                          :: setLabFrame
    procedure                          :: isInCMframe
  end type emissionENDF

abstract interface

  subroutine getAngleEnergy(self,angle,energy,rand )
    !! Interface for a subroutine of emissionsENDF class that returns angle and energy of emitted
    !! secondary neutron from a reaction given by class MT number.
    import :: defReal, &
              emissionENDF, &
              RNG
    class(emissionENDF), intent(in)   :: self
    real(defReal), intent(inout)      :: angle
    real(defReal), intent(inout)      :: energy
    class(RNG), intent(inout)         :: rand
  end subroutine

  subroutine getNumber(self,number)
    !! Interface for a subroutine of emissionsENDF class that returns average number of secondary
    !! neutrons emitted from the reactio ngiven by the class MT number. (0 for absorbtion)
    import :: defReal, &
              emissionENDF
    class(emissionENDF), intent(in) :: self
    real(defReal), intent(inout)     :: number
  end subroutine

end interface
contains

  subroutine setMT(self,MT)
    !! Subrutine that sets MT number of reaction the emissionENDF object is associated with
    class(emissionENDF), intent(inout) :: self
    integer(shortInt),intent(in)       :: MT
      self % MT = MT
  end subroutine setMT


  subroutine setLabFrame(self)
    class(emissionENDF), intent(inout) :: self

    self % cmFrame = .false.

  end subroutine setLabFrame


  function isInCMframe(self)
    class(emissionENDF), intent(in) :: self
    logical(defBool)                :: isInCMframe

    isInCMframe = self % cmFrame
  end function isInCMframe


end module emissionENDF_class
