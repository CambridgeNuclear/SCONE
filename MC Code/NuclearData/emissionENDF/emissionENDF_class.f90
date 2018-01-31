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
    procedure(sampleAngleEnergy),deferred :: sampleAngleEnergy
    procedure(releaseAt),deferred         :: releaseAt
    procedure                             :: setMT
    procedure                             :: setLabFrame
    procedure                             :: isInCMframe
  end type emissionENDF

  abstract interface

    subroutine sampleAngleEnergy(self,angle,E_out,E_in,rand )
      !! Interface for a subroutine of emissionsENDF class that returns angle and energy of emitted
      !! secondary neutron from a reaction given by class MT number.
      import :: defReal, &
                emissionENDF, &
                RNG
      class(emissionENDF), intent(in)   :: self
      real(defReal), intent(out)        :: angle
      real(defReal), intent(out)        :: E_out
      real(defReal), intent(in)         :: E_in
      class(RNG), intent(inout)         :: rand

    end subroutine sampleAngleEnergy

    function releaseAt(self,E_in) result(number)
      !! Interface for a subroutine of emissionsENDF class that returns average number of secondary
      !! neutrons emitted from the reactio ngiven by the class MT number. (0 for absorbtion)
      import :: defReal, &
                emissionENDF
      class(emissionENDF), intent(in)  :: self
      real(defReal), intent(in)        :: E_in
      real(defReal)                    :: number
    end function releaseAt

  end interface

  type,public :: emissionENDF_ptr
   !! Pointer Wrapper for emissionENDF class so arrays of pointers can be created
      private
      class(emissionENDF), pointer :: ptr => null()
    contains
      generic   :: assignment(=)     => assignPointer_ptr, assignPointer
      procedure :: sampleAngleEnergy => sampleAngleEnergy_ptr
      procedure :: releaseAt         => releaseAt_ptr
      procedure :: setMT             => setMT_ptr
      procedure :: setLabFrame       => setLabFrame_ptr
      procedure :: isInCMframe       => isInCMframe_ptr

      procedure, private :: assignPointer_ptr
      procedure, private :: assignPointer

  end type emissionENDF_ptr

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

!**************************************************************************************************!
! Pointer Wrapper Procedures
!
!**************************************************************************************************!


  subroutine setMT_ptr(self,MT)
    class(emissionENDF_ptr), intent(inout) :: self
    integer(shortInt),intent(in)           :: MT

    call self % ptr %setMT(MT)

  end subroutine setMT_ptr


  subroutine setLabFrame_ptr(self)
    class(emissionENDF_ptr), intent(inout) :: self

    call self % ptr % setLabFrame()

  end subroutine setLabFrame_ptr


  function isInCMframe_ptr(self)
    class(emissionENDF_ptr), intent(in) :: self
    logical(defBool)                    :: isInCMframe_ptr

    isInCMframe_ptr = self % ptr % isInCMframe()

  end function isInCMframe_ptr


  subroutine sampleAngleEnergy_ptr(self,angle,E_out,E_in,rand )
      class(emissionENDF_ptr), intent(in)   :: self
      real(defReal), intent(out)            :: angle
      real(defReal), intent(out)            :: E_out
      real(defReal), intent(in)             :: E_in
      class(RNG), intent(inout)             :: rand

      call self % ptr % sampleAngleEnergy(angle,E_out,E_in,rand)

  end subroutine sampleAngleEnergy_ptr

  function releaseAt_ptr(self,E_in) result(number)
    class(emissionENDF_ptr), intent(in) :: self
    real(defReal), intent(in)           :: E_in
    real(defReal)                       :: number

    number = self % ptr % releaseAt(E_in)

  end function releaseAt_ptr


  subroutine assignPointer_ptr(LHS,RHS)
    class(emissionENDF_ptr),intent(out) :: LHS
    type(emissionENDF_ptr),intent(in)   :: RHS

    if(associated(LHS % ptr)) deallocate(LHS % ptr)
    LHS % ptr => RHS % ptr

  end subroutine

  subroutine assignPointer(LHS,RHS)
    class(emissionENDF_ptr),intent(out)    :: LHS
    class(emissionENDF),pointer,intent(in) :: RHS

    if(associated(LHS % ptr)) deallocate(LHS % ptr)
    LHS % ptr => RHS
  end subroutine assignPointer




end module emissionENDF_class
