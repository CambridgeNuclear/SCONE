module correlatedReactionCE_inter

  use numPrecision
  use universalVariables
  use RNG_class,            only : RNG
  use reactionHandle_inter, only : reactionHandle

  implicit none
  private

  !!
  !! Public pointer cast
  !!
  public :: correlatedReactionCE_CptrCast

  !!
  !! Outgoing particle sample data
  !!
  !! Helper type to group data for a single outgoing particle sample
  !!
  !! Public members:
  !!   mu    -> cosine of polar deflection angle in [0;1]
  !!   phi   -> azimuthal deflection angle in [0;2*pi]
  !!   E_out -> outgoing energy of the particle [MeV]
  !!
  type, public :: outSample
    real(defReal) :: mu    = ZERO
    real(defReal) :: phi   = ZERO
    real(defReal) :: E_out = ZERO
  end type

  !!
  !! Multiple samples of outgoing particles
  !!
  !! Helper type to return an unknown number of correlated particles samples
  !! Uses global parameter MAX_OUTGOING_PARTICLES to know at compile-time how
  !! much space to reserve for particles data. It needs to avoid memory allocation
  !! on run-time for efficiency
  !!
  !! Public members:
  !!   N       -> number of samples generated N is integer in [0, MAX_OUTGOING_PARTICLES]
  !!   samples -> array of outSample type. Each enetry has data for a single outgoing sample
  !!
  type, public :: outgoingSamples
    integer(shortInt)                                 :: N = 0
    type(outSample),dimension(MAX_OUTGOING_PARTICLES) :: samples
  end type


  !!
  !! Reaction that produces secendary reaction that are correleated with each other
  !!
  !! This means that thare is a significant correlation between emitted particles.
  !! As the result all of outgoing particles need to be sampled at once
  !!
  !! Interface:
  !!   inCMframe -> returns true if reaction is in Centre-Of-Mass frame
  !!   sampleOut -> returns "outgoingSamples" with samples of outgoing particles
  !!
  type, public, abstract, extends(reactionHandle) :: correlatedReactionCE
    private
  contains
    procedure(inCMframe), deferred :: inCMframe
    procedure(sampleOut), deferred :: sampleOut
  end type correlatedReactionCE

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
      import :: defBool, correlatedReactionCE
      class(correlatedReactionCE), intent(in) :: self
      logical(defBool)                        :: isIt
    end function inCMframe

    !!
    !! Sample all outgoing particles at once
    !!
    !! Writes data into outgoingSamples type
    !! Sets number of outgoing particles N and writes data for particles 1-N
    !!
    !! Args:
    !!   particles [out] -> outgoingSamples type with space for data for outgoing particles
    !!   E_in [in]       -> incident particle energy [MeV]
    !!   rand [inout] -> Random number generator for samples
    !!
    !! Errors:
    !!   If E_in is out of bounds or invalid (e.g.-ve) fatalError is returned
    !!   If number of particles sampled N is larger than MAX_OUTGOING_PARTICLES
    !!   warning shall be produced (fatalError for now)
    !!
    !! TODO: Update documentation once warnings are implemented
    !!
    subroutine sampleOut(self, particles, E_in, rand)
      import :: outgoingSamples, defReal, correlatedReactionCE, RNG
      class(correlatedReactionCE), intent(in) :: self
      type(outgoingSamples), intent(out)      :: particles
      real(defReal), intent(in)               :: E_in
      class(RNG), intent(inout)               :: rand
    end subroutine sampleOut
  end interface

contains

  !!
  !! Cast reactionHandle pointer to correlatedReactionCE pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class reactionHandle
  !!
  !! Result:
  !!   Null is source is not of correlatedReactionCE class
  !!   Target points to source if source is correlatedReactionCE class
  !!
  pure function correlatedReactionCE_CptrCast(source) result(ptr)
    class(reactionHandle), pointer, intent(in) :: source
    class(correlatedReactionCE), pointer       :: ptr

    select type(source)
      class is(correlatedReactionCE)
        ptr => source

      class default
        ptr => null()
    end select

  end function correlatedReactionCE_CptrCast
    
end module correlatedReactionCE_inter
