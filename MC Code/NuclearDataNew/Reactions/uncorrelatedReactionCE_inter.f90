module uncorrelatedReactionCE_inter

  use numPrecision
  use RNG_class,            only : RNG
  use reactionHandle_inter, only : reactionHandle

  implicit none
  private


  !!
  !! Reaction that produces secendary reaction that are uncorreleated with each other
  !!
  !! This means that all particles produced are sampled from the same distribution.
  !! They are Indipendent and Identically Distributed.
  !!
  !! Interface:
  !!   inCMframe       -> returns true if reaction is in Centre-Of-Mass frame
  !!   release         -> return number (defReal) of particles producerd on average
  !!   releasePrompt   -> return number (defReal) of particles produced on average without delay
  !!   releaseDelayed  -> return number (defReal) of particles produced on average with delay
  !!   sampleDelayRate -> sample and return delay rate [1/s] assuming exponential distribution
  !!     for the decay.
  !!   sampleOutgoing  -> sample energy and the deflection of outgoing particle
  !!
  type, public, abstract, extends(reactionHandle) :: uncorrelatedReactionCE
    private
  contains

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
    pure function release(self, E) result(N)
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
    pure function releasePrompt(self, E) result(N)
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
    pure function releaseDelayed(self, E) result(N)
      import :: defReal, uncorrelatedReactionCE
      class(uncorrelatedReactionCE), intent(in) :: self
      real(defReal), intent(in)                 :: E
      real(defReal)                             :: N
    end function releaseDelayed

    !!
    !! Sample the delay rate for the delayed particle
    !!
    !! Assumption is made that random delay rate and particle distribution is independent.
    !! The delay rate assumes that the random neutron delay has exponential distribution
    !! (like radioactive decay). Units are [1/s]
    !!
    !! Args:
    !!   E [in]       -> Incident particle energy [MeV]
    !!   rand [inout] -> Random number generator for samples
    !!


  end interface


contains


end module uncorrelatedReactionCE_inter
