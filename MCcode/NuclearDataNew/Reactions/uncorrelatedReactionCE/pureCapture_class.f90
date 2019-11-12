module pureCapture_class

  use numPrecision
  use endfConstants
  use genericProcedures,            only : fatalError, numToChar
  use RNG_class,                    only : RNG
  use dataDeck_inter,               only : dataDeck
  use aceCard_class,                only : aceCard

  ! Interfaces
  use reactionHandle_inter,         only : reactionHandle
  use uncorrelatedReactionCE_inter, only : uncorrelatedReactionCE

  implicit none
  private

  !!
  !! Public pointer cast
  !!
  public :: pureCapture_TptrCast

  !!
  !! Pure capture reaction that produces so 2nd-ary particles
  !!
  !! To be used as a placeholder for kinematics for reactions without any 2nd-ary particles
  !!
  !! Interface:
  !!   uncorrelatedReactionCE interface
  !!
  type, public, extends(uncorrelatedReactionCE) :: pureCapture
    private

  contains
    !! Superclass interface
    procedure :: init
    procedure :: kill
    procedure :: inCMframe
    procedure :: release
    procedure :: releasePrompt
    procedure :: releaseDelayed
    procedure :: sampleOut
    procedure :: probOf
  end type pureCapture

contains

  !!
  !! Initialsie
  !!
  !! See reactionHandle for details
  !!
  !! Can be built for any reaction
  !!
  !! Errors:
  !!   None
  subroutine init(self, data, MT)
    class(pureCapture), intent(inout) :: self
    class(dataDeck), intent(inout)       :: data
    integer(shortInt), intent(in)        :: MT

  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  !! See reactionHandle for details
  !!
  elemental subroutine kill(self)
    class(pureCapture), intent(inout) :: self

  end subroutine kill

  !!
  !! Returns true if reaction is in Centre-Of-Mass frame
  !!
  !! See uncorrelatedReactionCE for details
  !!
  pure function inCMframe(self) result(isIt)
    class(pureCapture), intent(in) :: self
    logical(defBool)                  :: isIt

     isIt = .false.

  end function inCMframe

  !!
  !! Returns number of particles produced on average by the reaction
  !!
  !! See uncorrelatedReactionCE for details
  !!
  pure function release(self, E) result(N)
    class(pureCapture), intent(in) :: self
    real(defReal), intent(in)      :: E
    real(defReal)                  :: N

    N = ZERO

  end function release

  !!
  !! Returns number of particles produced on average instantly by the reaction
  !!
  !! See uncorrelatedReactionCE for details
  !!
  pure function releasePrompt(self, E) result(N)
    class(pureCapture), intent(in) :: self
    real(defReal), intent(in)      :: E
    real(defReal)                  :: N

    N = ZERO

  end function releasePrompt

  !!
  !! Returns number of particles produced on average by the reaction with delay
  !!
  !! See uncorrelatedReactionCE for details
  !!
  pure function releaseDelayed(self, E) result(N)
    class(pureCapture), intent(in) :: self
    real(defReal), intent(in)                :: E
    real(defReal)                            :: N

    N = ZERO

  end function releaseDelayed

  !!
  !! Sample outgoing particle
  !!
  !! See uncorrelatedReactionCE for details
  !!
  subroutine sampleOut(self, mu, phi, E_out, E_in, rand, lambda)
    class(pureCapture), intent(in)         :: self
    real(defReal), intent(out)             :: mu
    real(defReal), intent(out)             :: phi
    real(defReal), intent(out)             :: E_out
    real(defReal), intent(in)              :: E_in
    class(RNG), intent(inout)              :: rand
    real(defReal), intent(out), optional   :: lambda

    E_out = E_in
    mu = ONE
    phi = ZERO

    ! Default to the prompt particle
    if(present(lambda)) lambda = huge(lambda)

  end subroutine sampleOut

  !!
  !! Return probability density of emission at given angle and energy
  !!
  !! See uncorrelatedReactionCE for details
  !!
  function probOf(self, mu, phi, E_out, E_in) result(prob)
    class(pureCapture), intent(in) :: self
    real(defReal), intent(in)         :: mu
    real(defReal), intent(in)         :: phi
    real(defReal), intent(in)         :: E_out
    real(defReal), intent(in)         :: E_in
    real(defReal)                     :: prob

    prob = ONE
    if(E_out /= E_in)                           prob =  ZERO
    if(mu    /= ONE )                           prob =  ZERO
    if(.not.(phi   == ZERO .or. phi == TWO_PI)) prob = ZERO

  end function probOf

  !!
  !! Cast reactionHandle pointer to pureCapture pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class reactionHandle
  !!
  !! Result:
  !!   Null is source is not of pureCapture type
  !!   Target points to source if source is pureCapture type
  !!
  pure function pureCapture_TptrCast(source) result(ptr)
    class(reactionHandle), pointer, intent(in) :: source
    type(pureCapture), pointer                 :: ptr

    select type(source)
      type is(pureCapture)
        ptr => source

      class default
        ptr => null()

    end select

  end function pureCapture_TptrCast



end module pureCapture_class
