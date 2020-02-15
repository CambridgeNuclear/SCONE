module uncorrelatedReactionCE_inter

  use numPrecision
  use RNG_class,            only : RNG
  use reactionHandle_inter, only : reactionHandle

  implicit none
  private

  !!
  !! Public pointer cast
  !!
  public :: uncorrelatedReactionCE_CptrCast

  !!
  !! Reaction that produces secendary reaction that are uncorreleated with each other
  !!
  !! This means that all particles produced are sampled from the same distribution.
  !! They are Independent and Identically Distributed.
  !!
  !! Interface:
  !!   inCMframe       -> returns true if reaction is in Centre-Of-Mass frame
  !!   release         -> return number (defReal) of particles producerd on average
  !!   releasePrompt   -> return number (defReal) of particles produced on average without delay
  !!   releaseDelayed  -> return number (defReal) of particles produced on average with delay
  !!   sampleOut       -> sample energy and the deflection of outgoing particle
  !!   probOf          -> probability of emission at a given angle and energy
  !!
  type, public, abstract, extends(reactionHandle) :: uncorrelatedReactionCE
    private
  contains
    procedure(inCMframe),deferred       :: inCMframe
    procedure(release),deferred         :: release
    procedure(releasePrompt),deferred   :: releasePrompt
    procedure(releaseDelayed),deferred  :: releaseDelayed
    procedure(sampleOut), deferred      :: sampleOut
    procedure(probOf),deferred          :: probOf
  end type uncorrelatedReactionCE


  abstract interface

    !!
    !! Returns true if reaction is in Centre-Of-Mass frame
    !!
    !! Result:
    !!   .true. if in Centre-Of-Mass; .false. if in LAB-frame
    !!
    !! Errors:
    !!   None
    !!
    pure function inCMframe(self) result(isIt)
      import :: defBool, uncorrelatedReactionCE
      class(uncorrelatedReactionCE), intent(in) :: self
      logical(defBool)                          :: isIt
    end function inCMframe

    !!
    !! Returns number of particles produced on average by the reaction
    !!
    !! Args:
    !!   E [in] -> Incident particle energy [MeV]
    !!
    !! Result:
    !!   Average number of particles for an incident energy E.
    !!
    !! Errors:
    !!   If E is invalid (e.g. -ve) or outside of bounds of data table N = 0.0 is returned.
    !!
    function release(self, E) result(N)
      import :: defReal, uncorrelatedReactionCE
      class(uncorrelatedReactionCE), intent(in) :: self
      real(defReal), intent(in)                 :: E
      real(defReal)                             :: N
    end function release

    !!
    !! Returns number of particles produced on average instantly by the reaction
    !!
    !! Args:
    !!   E [in] -> Incident particle energy [MeV]
    !!
    !! Result:
    !!   Average number of prompt particles for an incident energy E.
    !!
    !! Errors:
    !!   If E is invalid (e.g. -ve) or outside of bounds of data table N = 0.0 is returned.
    !!
    function releasePrompt(self, E) result(N)
      import :: defReal, uncorrelatedReactionCE
      class(uncorrelatedReactionCE), intent(in) :: self
      real(defReal), intent(in)                 :: E
      real(defReal)                             :: N
    end function releasePrompt

    !!
    !! Returns number of particles produced on average by the reaction with delay
    !!
    !! Args:
    !!   E [in] -> Incident particle energy [MeV]
    !!
    !! Result:
    !!   Average number of delayed particles for an incident energy E.
    !!
    !! Errors:
    !!   If E is invalid (e.g. -ve) or outside of bounds of data table N = 0.0 is returned.
    !!
    function releaseDelayed(self, E) result(N)
      import :: defReal, uncorrelatedReactionCE
      class(uncorrelatedReactionCE), intent(in) :: self
      real(defReal), intent(in)                 :: E
      real(defReal)                             :: N
    end function releaseDelayed

    !!
    !! Sample outgoing particle
    !!
    !! Returns deflection from the incident direction and energy of an outgoing particle
    !! Uses either LAB or Centre-of-Mass reference frame.
    !!
    !! Can sample from both the prompt and the delayed particles.
    !! The delay rate assumes that the random neutron delay has exponential distribution
    !! (like radioactive decay). Units are [1/s].
    !! Delay rate of prompt particles is huge(defReal).
    !!
    !!
    !! Args:
    !!   mu  [out]    -> cosine of polar deflection angle in [0;1]
    !!   phi [out]    -> azimuthal deflection angle in [0;2*pi]
    !!   E_out [out]  -> outgoing particle energy [MeV]
    !!   E_in [in]    -> incident particle energy [MeV]
    !!   rand [inout] -> Random number generator for samples
    !!   lambda [out] -> Optional. Delay rate of the emission [1/s]
    !!
    !! Errors:
    !!   If E_in is out of bounds or invalid (e.g.-ve) fatalError is returned
    !!
    subroutine sampleOut(self, mu, phi, E_out, E_in, rand, lambda)
      import :: defReal, uncorrelatedReactionCE, RNG
      class(uncorrelatedReactionCE), intent(in) :: self
      real(defReal), intent(out)                :: mu
      real(defReal), intent(out)                :: phi
      real(defReal), intent(out)                :: E_out
      real(defReal), intent(in)                 :: E_in
      class(RNG), intent(inout)                 :: rand
      real(defReal), intent(out),optional       :: lambda
    end subroutine sampleOut

    !!
    !! Return probability density of emission at given angle and energy
    !!
    !! Probably not going to be used in calculation schemes very much but
    !! is very convinent for testing that the build resulted in correct
    !! probability distrbutiuon being loded
    !!
    !! Args:
    !!   mu    [in] -> outgoing cosine of polar deflection; mu in [0;1]
    !!   phi   [in] -> outgoing angle of azimuthal deflection; phi in [0,2*pi]
    !!   E_out [in] -> outgoing energy [MeV]
    !!   E_in  [in] -> incident energy [MeV]
    !!
    !! Result:
    !!   Value of probability distribution for the given arguments (must be +ve)
    !!
    !! Errors:
    !!   For input arguments outside their range returns prob = 0.0
    !!
    function probOf(self, mu, phi, E_out, E_in) result(prob)
      import :: defReal, uncorrelatedReactionCE
      class(uncorrelatedReactionCE), intent(in) :: self
      real(defReal), intent(in)                 :: mu
      real(defReal), intent(in)                 :: phi
      real(defReal), intent(in)                 :: E_out
      real(defReal), intent(in)                 :: E_in
      real(defReal)                             :: prob
    end function probOf
  end interface


contains

  !!
  !! Cast reactionHandle pointer to uncorrelatedReactionCE pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class reactionHandle
  !!
  !! Result:
  !!   Null is source is not of uncorrelatedReactionCE class
  !!   Target points to source if source is uncorrelatedReactionCE class
  !!
  pure function uncorrelatedReactionCE_CptrCast(source) result(ptr)
    class(reactionHandle), pointer, intent(in) :: source
    class(uncorrelatedReactionCE), pointer     :: ptr

    select type(source)
      class is(uncorrelatedReactionCE)
        ptr => source

      class default
        ptr => null()
    end select

  end function uncorrelatedReactionCE_CptrCast

end module uncorrelatedReactionCE_inter
