module elasticNeutronScatter_class

  use numPrecision
  use genericProcedures,            only : fatalError
  use RNG_class,                    only : RNG
  use dataDeck_inter,               only : dataDeck
  use aceCard_class,                only : aceCard
  use uncorrelatedReactionCE_inter, only : uncorrelatedReactionCE


  implicit none
  private

  !!
  !! Reaction type for Neutron Elastic Scattering
  !!
  !! Implements standard elastic neutron stattering
  !!
  !! Interface:
  !!   uncorrelatedReactionCE interface
  !!   buildFromACE -> initialise object from ACE dataCard
  type, public, extends(uncorrelatedReactionCE) :: elasticNeutronScatter
    private

  contains
    !! Superclass interface
    procedure :: init
    procedure :: kill
    procedure :: inCMframe
    procedure :: release
    procedure :: releasePrompt
    procedure :: releaseDelayed
    procedure :: sampleDelayRate
    procedure :: sampleOut

    !! Instance procedures
    procedure :: buildFromACE

  end type elasticNeutronScatter

contains

  !!
  !! Initialsie
  !!
  !! See reactionHandle for details
  !!
  !! Errors:
  !!   fatalError if MT is not neutron elastic scattering
  !!
  subroutine init(self, data, MT)
    class(elasticNeutronScatter), intent(inout) :: self
    class(dataDeck), intent(inout)              :: data
    integer(shortInt), intent(in)               :: MT
    character(100), parameter :: Here = 'init (elasticNeutronScatter_class.f90)'

    select type(data)
      type is (aceCard)
        call self % buildFromACE(data)

      class default
        call fatalError(Here,'Elastic neutron scattering cannot be build from '//data % myType())
    end select

  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  !! See reactionHandle for details
  !!
  elemental subroutine kill(self)
    class(elasticNeutronScatter), intent(inout) :: self

    ! Nothing to do yet

  end subroutine kill

  !!
  !! Returns true if reaction is in Centre-Of-Mass frame
  !!
  !! See uncorrelatedReactionCE for details
  !!
  pure function inCMframe(self) result(isIt)
    class(elasticNeutronScatter), intent(in) :: self
    logical(defBool)                         :: isIt

     isIt = .true.

  end function inCMframe

  !!
  !! Returns number of particles produced on average by the reaction
  !!
  !! See uncorrelatedReactionCE for details
  !!
  pure function release(self, E) result(N)
    class(elasticNeutronScatter), intent(in) :: self
    real(defReal), intent(in)                :: E
    real(defReal)                            :: N

    N = ONE

  end function release

  !!
  !! Returns number of particles produced on average instantly by the reaction
  !!
  !! See uncorrelatedReactionCE for details
  !!
  pure function releasePrompt(self, E) result(N)
    class(elasticNeutronScatter), intent(in) :: self
    real(defReal), intent(in)                :: E
    real(defReal)                            :: N

    N = ONE

  end function releasePrompt
    
  !!
  !! Returns number of particles produced on average by the reaction with delay
  !!
  !! See uncorrelatedReactionCE for details
  !!
  pure function releaseDelayed(self, E) result(N)
    class(elasticNeutronScatter), intent(in) :: self
    real(defReal), intent(in)                :: E
    real(defReal)                            :: N

    N = ZERO

  end function releaseDelayed

  !!
  !! Sample the delay rate for the delayed particle
  !!
  !! See uncorrelatedReactionCE for details
  !!
  function sampleDelayRate(self, E, rand) result(lambda)
    class(elasticNeutronScatter), intent(in) :: self
    real(defReal), intent(in)                :: E
    class(RNG), intent(inout)                :: rand
    real(defReal)                            :: lambda

    lambda = ZERO

  end function sampleDelayRate

  !!
  !! Sample outgoing particle
  !!
  !! See uncorrelatedReactionCE for details
  !!
  subroutine sampleOut(self, mu, phi, E_out, E_in, rand)
    class(elasticNeutronScatter), intent(in) :: self
    real(defReal), intent(out)               :: mu
    real(defReal), intent(out)               :: phi
    real(defReal), intent(out)               :: E_out
    real(defReal), intent(in)                :: E_in
    class(RNG), intent(inout)                :: rand

    mu  = ONE
    phi = ZERO
    E_out = E_in

  end subroutine sampleOut

  !!
  !! Build elasticNeutronScatter from ACE dataCard
  !!
  subroutine buildFromACE(self, ACE)
    class(elasticNeutronScatter), intent(inout) :: self
    type(aceCard), intent(inout)                :: ACE

    ! Write it

  end subroutine buildFromACE

end module elasticNeutronScatter_class
