module reactionMG_inter

  use numPrecision
  use RNG_class,            only : RNG
  use reactionHandle_inter, only : reactionHandle

  implicit none
  private


  !!
  !! Public pointer cast
  !!
  public :: reactionMG_ptrCast

  !!
  !! Reaction that produces secendary particles in MG
  !!
  !! It is assumed that all outgoing particles are sampled from the same distribution.
  !! They are Independent and Identically Distributed.
  !! Furthermore reactions are assumed to take place in LAB reference frame
  !! In fact Centre-Of-Mass frame has little use for MG reactions as transition
  !! between CoM and LAB frames cannoct be done without intermediate conversion to CE.
  !!
  !! Interface:
  !!   release         -> return number (defReal) of particles producerd on average
  !!   releasePrompt   -> return number (defReal) of particles produced on average without delay
  !!   releaseDelayed  -> return number (defReal) of particles produced on average with delay
  !!   sampleDelayRate -> sample and return delay rate [1/s] assuming exponential distribution
  !!     for the decay.
  !!   sampleOut       -> sample energy and the deflection of outgoing particle
  !!
  type, public, abstract, extends(reactionHandle) :: reactionMG
    private
  contains
    procedure(release),deferred         :: release
    procedure(releasePrompt),deferred   :: releasePrompt
    procedure(releaseDelayed),deferred  :: releaseDelayed
    procedure(sampleDelayRate),deferred :: sampleDelayRate
    procedure(sampleOut), deferred      :: sampleOut
  end type reactionMG

  abstract interface

    !!
    !! Returns number of particles produced on average by the reaction
    !!
    !! Args:
    !!   G [in] -> Incident particle energy group
    !!
    !! Result:
    !!   Average number of particles for an incident energy group G.
    !!
    !! Errors:
    !!   If G is invalid (e.g. -ve) or outside of bounds of data table N = 0.0 is returned.
    !!
    pure function release(self, G) result(N)
      import :: defReal, shortInt, reactionMG
      class(reactionMG), intent(in) :: self
      integer(shortInt),intent(in)  :: G
      real(defReal)                 :: N
    end function release

    !!
    !! Returns number of particles produced instantly by the reaction
    !!
    !! Args:
    !!   G [in] -> Incident particle energy group
    !!
    !! Result:
    !!   Average number of prompt particles for an incident energy group G.
    !!
    !! Errors:
    !!   If G is invalid (e.g. -ve) or outside of bounds of data table N = 0.0 is returned.
    !!
    pure function releasePrompt(self, G) result(N)
      import :: defReal, shortInt, reactionMG
      class(reactionMG), intent(in) :: self
      integer(shortInt),intent(in)  :: G
      real(defReal)                 :: N
    end function releasePrompt

    !!
    !! Returns number of particles produced instantly by the reaction
    !!
    !! Args:
    !!   G [in] -> Incident particle energy group
    !!
    !! Result:
    !!   Average number of prompt particles for an incident energy group G.
    !!
    !! Errors:
    !!   If G is invalid (e.g. -ve) or outside of bounds of data table N = 0.0 is returned.
    !!
    pure function releaseDelayed(self, G) result(N)
      import :: defReal, shortInt, reactionMG
      class(reactionMG), intent(in) :: self
      integer(shortInt),intent(in)  :: G
      real(defReal)                 :: N
    end function releaseDelayed

    !!
    !! Sample the delay rate for the delayed particle
    !!
    !! Assumption is made that random delay rate and particle distribution is independent.
    !! The delay rate assumes that the random neutron delay has exponential distribution
    !! (like radioactive decay). Units are [1/s]
    !!
    !! Args:
    !!   G [in]       -> Incident particle energy group
    !!   rand [inout] -> Random number generator for samples
    !!
    !! Result:
    !!   Time constant (lambda) of exponentialy distributed neutron emission delay in [1/s]
    !!
    !! Errors:
    !!   If the reaction has no delayed emissions lambda = 0.0 is returned.
    !!   If the reaction does not produce delayed emissions at a given energy,
    !!   energy group G is invalid (e.g. -ve) or outside bounds of stored data lambda = 0.0
    !!   is returned
    !!
    function sampleDelayRate(self, G, rand) result(lambda)
      import :: defReal, ,shortInt, reactionMG, RNG
      class(reactionMG), intent(in) :: self
      integer(shortInt), intent(in) :: G
      class(RNG), intent(inout)     :: rand
    end function sampleDelayRate

    !!
    !! Sample outgoing particle
    !!
    !! Returns deflection from the incident direction and energy group of an outgoing particle
    !!
    !! Args:
    !!   mu  [out]    -> cosine of polar deflection angle in [0;1]
    !!   phi [out]    -> azimuthal deflection angle in [0;2*pi]
    !!   G_out [out]  -> outgoing particle energy group
    !!   G_in [in]    -> incident particle energy group
    !!   rand [inout] -> Random number generator for samples
    !!
    !! Errors:
    !!   If G_in is out of bounds or invalid (e.g.-ve) fatalError is returned
    !!
    subroutine sampleOut(self, mu, phi, G_out, G_in, rand)
      import :: defReal, shortInt, reactionMG, RNG
      class(reactionMG), intent(in)  :: self
      real(defReal), intent(out)     :: mu
      real(defReal), intent(out)     :: phi
      integer(shortInt), intent(out) :: G_out
      integer(shortInt), intent(in)  :: G_in
      class(RNG), intent(inout)      :: rand
    end subroutine sampleOut
  end interface

contains

  !!
  !! Cast reactionHandle pointer to reactionMG pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class reactionHandle
  !!
  !! Result:
  !!   Null is source is not of reactionMG class
  !!   Target points to source if source is reactionMG class
  !!
  !! NOTE:
  !!   If target is a unique reference to an object. Memory leak will be caused
  !!
  pure function reactionMG_ptrCast(source) result(ptr)
    class(reactionHandle), pointer, intent(in) :: source
    class(reactionMG), pointer                 :: ptr

    select type(source)
      class is(reactionMG)
        ptr => source

      class default
        ptr => null()
    end select

  end function reactionMG_ptrCast

end module reactionMG_inter
