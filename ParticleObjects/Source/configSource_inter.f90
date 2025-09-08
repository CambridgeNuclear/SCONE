module configSource_inter

  use numPrecision
  use RNG_class,      only : RNG
  use particle_class, only : particleState
  use source_inter,   only : source,  kill_super => kill
  implicit none

  !!
  !! Extendable scource class procedures
  !!
  public :: kill

  !!
  !! Configurable source
  !!
  !! Generates a single particle sample by calling a number of subroutines
  !! related to each component
  !!
  !! A single sample of a particle state is created by calling `sample****` procedures in order
  !! given in `sampleParticle` function. Note that `sampleEnergyAngle` is called AFTER
  !! `sampleEnergy`.
  !!
  !! Interface:
  !!   source_inter Interface
  !!   sampleType        -> sets the particle type
  !!   samplePosition    -> samples the particle's position in the geometry
  !!   sampleEnergy      -> samples the particle's energy
  !!   sampleEnergyAngle -> samples the particle's energy and angle from corresponding distr.
  !!
  type, public, abstract, extends(source) :: configSource

  contains
    procedure                              :: sampleParticle
    procedure(sampleType), deferred        :: sampleType
    procedure(samplePosition), deferred    :: samplePosition
    procedure(sampleEnergy), deferred      :: sampleEnergy
    procedure(sampleEnergyAngle), deferred :: sampleEnergyAngle
  end type configSource

  abstract interface

    !!
    !! Sample Type of a particle
    !!
    !! Sets 'Type' in the particleState p (e.g. P_NEUTRON)
    !!
    !! Inputs:
    !!   p [inout] -> particleState to be given a type
    !!   rand [in] -> random number generator
    !!
    subroutine sampleType(self, p, rand)
      import :: configSource, &
                particleState, &
                RNG
      class(configSource), intent(inout)  :: self
      class(particleState), intent(inout) :: p
      class(RNG), intent(inout)           :: rand
    end subroutine sampleType

    !!
    !! Sample particle position
    !!
    !! Sets position of the particle p
    !!
    !! Inputs:
    !!   p [inout] -> particleState to be given a position
    !!   rand [in] -> random number generator
    !!
    subroutine samplePosition(self, p, rand)
      import :: configSource, &
                particleState, &
                RNG
      class(configSource), intent(inout)  :: self
      class(particleState), intent(inout) :: p
      class(RNG), intent(inout)           :: rand
    end subroutine samplePosition

    !!
    !! Sample particle Energy/Group
    !!
    !! Sets energy of a particle to a CE value or a MG index
    !! Also sets 'isMG' flag to .true. or .false.
    !!
    !! Inputs:
    !!   p [inout] -> particleState to be given a position
    !!   rand [in] -> random number generator
    !!
    subroutine sampleEnergy(self, p, rand)
      import :: configSource, &
                particleState, &
                RNG
      class(configSource), intent(inout)  :: self
      class(particleState), intent(inout) :: p
      class(RNG), intent(inout)           :: rand
    end subroutine sampleEnergy

    !!
    !! Sample particle Energy/Group and angle
    !!
    !! Sets direction of a particle together with its energy.
    !! Sampling of energy is optional if Angle & Energy are uncorrelated
    !! Is called after `sampleEnergy`, to overwrite value provided by that subroutine
    !!
    !! Inputs:
    !!   p [inout] -> particleState to be given a position
    !!   rand [in] -> random number generator
    !!
    subroutine sampleEnergyAngle(self, p, rand)
      import :: configSource, &
                particleState, &
                RNG
      class(configSource), intent(inout)  :: self
      class(particleState), intent(inout) :: p
      class(RNG), intent(inout)           :: rand
    end subroutine sampleEnergyAngle

  end interface

contains

  !!
  !! Sample particle's phase space co-ordinates
  !!
  !! See source_inter for details
  !!
  function sampleParticle(self, rand) result(p)
    class(configSource), intent(inout) :: self
    class(RNG), intent(inout)          :: rand
    type(particleState)                :: p

    call self % sampleType(p, rand)
    call self % samplePosition(p, rand)
    call self % sampleEnergyAngle(p, rand)
    call self % sampleEnergy(p, rand)
    p % time = ZERO
    p % wgt  = ONE

  end function sampleParticle

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(configSource), intent(inout) :: self

    call kill_super(self)

  end subroutine kill

end module configSource_inter
