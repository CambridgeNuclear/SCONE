module elasticNeutronScatter_class

  use numPrecision
  use endfConstants
  use genericProcedures,            only : fatalError, numToChar
  use RNG_class,                    only : RNG
  use dataDeck_inter,               only : dataDeck
  use aceCard_class,                only : aceCard
  use tabularAngle_class,           only : tabularAngle
  use reactionHandle_inter,         only : reactionHandle
  use uncorrelatedReactionCE_inter, only : uncorrelatedReactionCE

  implicit none
  private

  !!
  !! Public pointer cast
  !!
  public :: elasticNeutronScatter_TptrCast

  !!
  !! Reaction type for Neutron Elastic Scattering
  !!
  !! Implements standard elastic neutron stattering
  !!
  !! Interface:
  !!   uncorrelatedReactionCE interface
  !!   buildFromACE -> initialise object from ACE dataCard
  !!
  type, public, extends(uncorrelatedReactionCE) :: elasticNeutronScatter
    private
    type(tabularAngle) :: angularData
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
    procedure :: probOf

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

    ! Catch MT numbers that are not elastic scattering
    if( MT /= N_N_ELASTIC) then
      call fatalError(Here,'Connot be build for reaction with MT:' //numToChar(MT)//&
                           ' which is not scattering')
    end if

    ! Select buld procedure approperiate for given dataDeck
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

    call self % angularData % kill()

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

    ! Set energy
    E_out = E_in

    ! Sample mu
    mu = self % angularData % sample(E_in, rand)

    ! Sample phi
    phi = rand % get() * TWO_PI

  end subroutine sampleOut

  !!
  !! Return probability density of emission at given angle and energy
  !!
  !! See uncorrelatedReactionCE for details
  !!
  function probOf(self, mu, phi, E_out, E_in) result(prob)
    class(elasticNeutronScatter), intent(in) :: self
    real(defReal), intent(in)                :: mu
    real(defReal), intent(in)                :: phi
    real(defReal), intent(in)                :: E_out
    real(defReal), intent(in)                :: E_in
    real(defReal)                            :: prob

    ! Check range and set mu prob
    if (abs(mu) <= ONE .and. E_out > ZERO) then
      prob = self % angularData % probabilityOf(mu, E_in)
    else
      prob = ZERO
    end if

    ! Apply E_out prob -> delta distribution
    if( E_out /= E_in) prob = ZERO

    ! Apply phi prob
    if (phi >= ZERO .and. phi <= TWO_PI) then
      prob = prob /TWO_PI
    else
      prob = ZERO
    end if

  end function probOf

  !!
  !! Build elasticNeutronScatter from ACE dataCard
  !!
  subroutine buildFromACE(self, ACE)
    class(elasticNeutronScatter), intent(inout) :: self
    type(aceCard), intent(inout)                :: ACE

    ! Seat read head of ACE to elastic Scattering data
    call ACE % setToAngleMT(N_N_ELASTIC)

    ! Initialise tabular angle
    call self % angularData % init(ACE, N_N_ELASTIC)

  end subroutine buildFromACE

  !!
  !! Cast reactionHandle pointer to elasticNeutronScatter pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class reactionHandle
  !!
  !! Result:
  !!   Null is source is not of elasticNeutronScatter type
  !!   Target points to source if source is elasticNeutronScatter type
  !!
  pure function elasticNeutronScatter_TptrCast(source) result(ptr)
    class(reactionHandle), pointer, intent(in) :: source
    type(elasticNeutronScatter), pointer       :: ptr

    select type(source)
      type is(elasticNeutronScatter)
        ptr => source

      class default
        ptr => null()
    end select

  end function elasticNeutronScatter_TptrCast


end module elasticNeutronScatter_class
