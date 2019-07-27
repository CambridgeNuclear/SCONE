module fissionCE_class

  use numPrecision
  use endfConstants
  use genericProcedures,            only : fatalError, numToChar
  use RNG_class,                    only : RNG
  use dataDeck_inter,               only : dataDeck
  use aceCard_class,                only : aceCard
  use reactionHandle_inter,         only : reactionHandle
  use releaseLawENDF_inter,         only : releaseLawENDF
  use energyLawENDF_inter,          only : energyLawENDF
  use uncorrelatedReactionCE_inter, only : uncorrelatedReactionCE

  ! Factiories
  use releaseLawENDFfactory_func,   only : new_releaseLawENDF
  use energyLawENDFfactory_func,    only : new_energyLawENDF


  implicit none
  private

  !!
  !! Public pointer cast
  !!
  public :: fissionCE_TptrCast

  !!
  !! Simple form of a fission reaction used in eigenvalue calculations
  !!
  !! Does not contain any information about prompt and delayed neutrons
  !!
  !! Private Members:
  !!   nuBar -> energy dependant average neutron release
  !!   eLaw -> energy distributions for outgoing particles
  !!
  !! Interface:
  !!   uncorrelatedReactionCE interface
  !!
  type, public, extends(uncorrelatedReactionCE) :: fissionCE
    private
    class(releaseLawENDF),allocatable :: nuBar
    class(energyLawENDF),allocatable  :: eLaw
  contains
    ! Superclass procedures
    procedure :: init
    procedure :: kill
    procedure :: inCMframe
    procedure :: release
    procedure :: releasePrompt
    procedure :: releaseDelayed
    procedure :: sampleDelayRate
    procedure :: sampleOut
    procedure :: probOf

    ! Type specific procedures
    procedure :: buildFromACE
  end type fissionCE

contains

  !!
  !! Initialsie
  !!
  !! See reactionHandle for details
  !!
  !! Errors:
  !!   fatalError if MT /= N_FISSION
  subroutine init(self, data, MT)
    class(fissioNCE), intent(inout) :: self
    class(dataDeck), intent(inout)  :: data
    integer(shortInt), intent(in)   :: MT
    character(100),parameter :: Here ='init (fissionCE_class.f90)'

    if( MT /= N_FISSION) then
      call fatalError(Here,'fissionCE suports only MT=18. Was given: '//numToChar(MT))
    end if

    ! Select buld procedure approperiate for given dataDeck
    select type(data)
      type is (aceCard)
        call self % buildFromACE(data)

      class default
        call fatalError(Here,'Fission CE cannot be build from '//data % myType())
    end select


  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(fissionCE), intent(inout) :: self

    if(allocated(self % nuBar)) then
      call self % nuBar % kill()
      deallocate(self % nuBar)
    end if

    if(allocated(self % eLaw)) then
      call self % eLaw % kill()
      deallocate(Self % eLaw)
    end if

  end subroutine kill

  !!
  !! Returns true if reaction is in Centre-Of-Mass frame
  !!
  !! See uncorrelatedReactionCE for details
  !!
  pure function inCMframe(self) result(isIt)
    class(fissionCE), intent(in) :: self
    logical(defBool)             :: isIt

     isIt = .false.

  end function inCMframe

  !!
  !! Returns number of particles produced on average by the reaction
  !!
  !! See uncorrelatedReactionCE for details
  !!
  function release(self, E) result(N)
    class(fissionCE), intent(in) :: self
    real(defReal), intent(in)    :: E
    real(defReal)                :: N

    N = self % nuBar % releaseAt(E)

  end function release
    
  !!
  !! Returns number of particles produced on average instantly by the reaction
  !!
  !! See uncorrelatedReactionCE for details
  !!
  function releasePrompt(self, E) result(N)
    class(fissionCE), intent(in) :: self
    real(defReal), intent(in)    :: E
    real(defReal)                :: N

    N = self % nuBar % releaseAt(E)

  end function releasePrompt

  !!
  !! Returns number of particles produced on average by the reaction with delay
  !!
  !! See uncorrelatedReactionCE for details
  !!
  pure function releaseDelayed(self, E) result(N)
    class(fissionCE), intent(in) :: self
    real(defReal), intent(in)    :: E
    real(defReal)                :: N

    N = ZERO

  end function releaseDelayed

  !!
  !! Sample the delay rate for the delayed particle
  !!
  !! See uncorrelatedReactionCE for details
  !!
  function sampleDelayRate(self, E, rand) result(lambda)
    class(fissionCE), intent(in) :: self
    real(defReal), intent(in)    :: E
    class(RNG), intent(inout)    :: rand
    real(defReal)                :: lambda

    lambda = ZERO

  end function sampleDelayRate

  !!
  !! Sample outgoing particle
  !!
  !! See uncorrelatedReactionCE for details
  !!
  subroutine sampleOut(self, mu, phi, E_out, E_in, rand)
    class(fissionCE), intent(in) :: self
    real(defReal), intent(out)   :: mu
    real(defReal), intent(out)   :: phi
    real(defReal), intent(out)   :: E_out
    real(defReal), intent(in)    :: E_in
    class(RNG), intent(inout)    :: rand

    ! Sample mu
    mu = TWO * rand % get() - ONE

    ! Sample Phi
    phi = TWO_PI * rand % get()

    ! Sample E_out
    E_out = self % eLaw % sample(E_in, rand)

  end subroutine sampleOut

  !!
  !! Return probability density of emission at given angle and energy
  !!
  !! See uncorrelatedReactionCE for details
  !!
  function probOf(self, mu, phi, E_out, E_in) result(prob)
    class(fissionCE), intent(in) :: self
    real(defReal), intent(in)                :: mu
    real(defReal), intent(in)                :: phi
    real(defReal), intent(in)                :: E_out
    real(defReal), intent(in)                :: E_in
    real(defReal)                            :: prob

    if(abs(mu) <= ONE .and. E_out > ZERO .and. phi <= TWO_PI .and. phi >= ZERO) then
      prob = self % eLaw % probabilityOf(E_out, E_in) / (TWO * TWO_PI)
    else
      prob = ZERO
    end if

  end function probOf

  !!
  !! Build fissionCE from ACE dataCard
  !!
  subroutine buildFromACE(self, ACE)
    class(fissionCE), intent(inout) :: self
    type(aceCard), intent(inout)    :: ACE

    ! Read Release data
    call new_releaseLawENDF(self % nuBar, ACE, N_FISSION)

    ! Read Energy data
    call new_energyLawENDF(self % eLaw, ACE, N_FISSION)

  end subroutine buildFromACE

  !!
  !! Cast reactionHandle pointer to fissionCE pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class reactionHandle
  !!
  !! Result:
  !!   Null is source is not of fissionCE type
  !!   Target points to source if source is fissionCE type
  !!
  pure function fissionCE_TptrCast(source) result(ptr)
    class(reactionHandle), pointer, intent(in) :: source
    type(fissionCE), pointer                   :: ptr

    select type(source)
      type is(fissionCE)
        ptr => source

      class default
        ptr => null()

    end select

  end function fissionCE_TptrCast

end module fissionCE_class
