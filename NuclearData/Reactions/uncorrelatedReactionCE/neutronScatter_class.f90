module neutronScatter_class

  use numPrecision
  use endfConstants
  use genericProcedures,            only : fatalError, numToChar
  use RNG_class,                    only : RNG
  use dataDeck_inter,               only : dataDeck
  use aceCard_class,                only : aceCard

  ! Interfaces
  use reactionHandle_inter,         only : reactionHandle
  use uncorrelatedReactionCE_inter, only : uncorrelatedReactionCE

  ! ENDF Data Interfaces
  use angleLawENDF_inter,            only : angleLawENDF
  use energyLawENDF_inter,           only : energyLawENDF
  use correlatedLawENDF_inter,       only : correlatedLawENDF

  ! Factories
  use angleLawENDFfactory_func,      only : new_angleLawENDF
  use energyLawENDFfactory_func,     only : new_energyLawENDF
  use correlatedLawENDFfactory_func, only : new_correlatedLawENDF


  implicit none
  private

  !!
  !! Public pointer cast
  !!
  public :: neutronScatter_TptrCast

  !!
  !! Reaction Type for any neutron scattering
  !!
  !! This includes any reaction that produces secondary neutrons but is not fission
  !! For example : inelastic scattering, (N,2N), (N,alphaN) etc.
  !!
  !! Number of secondary emissions is independent of incident energy.
  !! No delayed emissions are supported
  !!
  !! Private Members:
  !!   cmFrame    -> true if reaction is in CM frame
  !!   correlated -> true if reaction is using correlated Energy-Angle Law
  !!   N_out      -> Average number of outgoing particles
  !!   muLaw      -> Angular law for outgoing neutron
  !!   eLaw       -> Energy law for outgoing neutron
  !!   corrLaw    -> Correlated low for outgoing neutron
  !!
  !! Interface:
  !!   uncorrelatedReactionCE interface
  !!   buildFromACE -> initialise object from ACE dataCard
  !!
  type, public, extends(uncorrelatedReactionCE) :: neutronScatter
    private
    ! Reference frame and emission type flags
    logical(defBool)   :: cmFrame    = .true.
    logical(defBool)   :: correlated = .false.
    real(defReal)      :: N_out      = ZERO

    ! Emission laws
    class(angleLawENDF),allocatable      :: muLaw
    class(energyLawENDF),allocatable     :: eLaw
    class(correlatedLawENDF),allocatable :: corrLaw
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

    !! Instance procedures
    procedure :: buildFromACE

  end type neutronScatter

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
    class(neutronScatter), intent(inout) :: self
    class(dataDeck), intent(inout)              :: data
    integer(shortInt), intent(in)               :: MT
    character(100), parameter :: Here = 'init (neutronScatter_class.f90)'

    ! Select buld procedure approperiate for given dataDeck
    select type(data)
      type is (aceCard)
        call self % buildFromACE(data, MT)

      class default
        call fatalError(Here,'Neutron scattering cannot be build from '//data % myType())
    end select

  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  !! See reactionHandle for details
  !!
  elemental subroutine kill(self)
    class(neutronScatter), intent(inout) :: self

    ! Reset constants
    self % cmFrame    = .true.
    self % correlated = .false.
    self % N_out      = ZERO

    ! Kill ENDF distributions
    if(allocated(self % muLaw))   call self % muLaw % kill()
    if(allocated(self % eLaw))    call self % eLaw % kill()
    if(allocated(self % corrLaw)) call self % corrLaw % kill()

    ! Deallocate ENDF distributions
    if(allocated(self % muLaw))   deallocate(self % muLaw)
    if(allocated(self % eLaw))    deallocate(self % eLaw)
    if(allocated(self % corrLaw)) deallocate(self % corrLaw)

  end subroutine kill

  !!
  !! Returns true if reaction is in Centre-Of-Mass frame
  !!
  !! See uncorrelatedReactionCE for details
  !!
  pure function inCMframe(self) result(isIt)
    class(neutronScatter), intent(in) :: self
    logical(defBool)                         :: isIt

     isIt = self % cmFrame

  end function inCMframe

  !!
  !! Returns number of particles produced on average by the reaction
  !!
  !! See uncorrelatedReactionCE for details
  !!
  pure function release(self, E) result(N)
    class(neutronScatter), intent(in) :: self
    real(defReal), intent(in)                :: E
    real(defReal)                            :: N

    N = self % N_out

  end function release

  !!
  !! Returns number of particles produced on average instantly by the reaction
  !!
  !! See uncorrelatedReactionCE for details
  !!
  pure function releasePrompt(self, E) result(N)
    class(neutronScatter), intent(in) :: self
    real(defReal), intent(in)                :: E
    real(defReal)                            :: N

    N = self % N_out

  end function releasePrompt

  !!
  !! Returns number of particles produced on average by the reaction with delay
  !!
  !! See uncorrelatedReactionCE for details
  !!
  pure function releaseDelayed(self, E) result(N)
    class(neutronScatter), intent(in) :: self
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
    class(neutronScatter), intent(in) :: self
    real(defReal), intent(out)               :: mu
    real(defReal), intent(out)               :: phi
    real(defReal), intent(out)               :: E_out
    real(defReal), intent(in)                :: E_in
    class(RNG), intent(inout)                :: rand
    real(defReal), intent(out), optional     :: lambda

    ! Sample energy an angle
    if( self % correlated) then
      call self % corrLaw % sample(mu, E_out, E_in, rand)

    else
      mu    = self % muLaw % sample(E_in,rand)
      E_out = self % eLaw  % sample(E_in,rand)

    end if

    ! Sample phi
    phi = rand % get() * TWO_PI

    ! Only prompt particles. Set delay
    if(present(lambda)) lambda = huge(lambda)

  end subroutine sampleOut

  !!
  !! Return probability density of emission at given angle and energy
  !!
  !! See uncorrelatedReactionCE for details
  !!
  function probOf(self, mu, phi, E_out, E_in) result(prob)
    class(neutronScatter), intent(in) :: self
    real(defReal), intent(in)         :: mu
    real(defReal), intent(in)         :: phi
    real(defReal), intent(in)         :: E_out
    real(defReal), intent(in)         :: E_in
    real(defReal)                     :: prob

    ! Check range
    if (abs(mu) <= ONE .and. E_out > ZERO) then
      if(self % correlated) then
        prob = self % corrLaw % probabilityOf(mu, E_out, E_in)
      else
        prob = self % muLaw % probabilityOf(mu, E_in)
        prob = prob * self % eLaw % probabilityOf(E_out, E_in)
      end if
    else
      prob = ZERO
    end if

    ! Apply phi prob
    if (phi >= ZERO .and. phi <= TWO_PI) then
      prob = prob /TWO_PI
    else
      prob = ZERO
    end if

  end function probOf

  !!
  !! Build neutronScatter from ACE dataCard
  !!
  !! Args:
  !!   ACE [inout] -> ACE Card class. Head set to any place
  !!   MT [in]     -> MT number of reaction requested
  !!
  !! Errors:
  !!   FatalError if MT number does not produce 2nd-ary neutrons
  !!   FatalError if MT number has energy dependant neutron yield
  !!
  subroutine buildFromACE(self, ACE, MT)
    class(neutronScatter), intent(inout) :: self
    type(aceCard), intent(inout)         :: ACE
    integer(shortInt), intent(in)        :: MT
    integer(shortInt)                    :: LOCB, TY
    character(100),parameter :: Here ='buildFromACE (neutronScatter_class.f90)'

    if( ACE % isCaptureMT(MT)) then
      call fatalError(Here, 'Requested reaction with MT: '//numToChar(MT)//&
                             ' does not produce 2nd-ary neutrons')
    end if

    ! Read LOCB for reaction under MT
    LOCB = ACE % LOCBforMT(MT)

    ! Read if data for reaction is in Centre-of-Mass frame
    self % cmFrame = ACE % isCMframe(MT)

    ! Read number of 2nd-ary particles
    TY = ACE % neutronReleaseMT(MT)
    if(TY == 19 .or. TY > 100) then
      call fatalError(Here,'Reaction with MT: '// numToChar(MT)//&
                           ' has energy dependent neutron yield. It is not supported')
    elseif(TY < 0) then
      call fatalError(Here,'-ve Neutron Release WTF?')

    end if

    self % N_out = real(TY, defReal)

    ! Build as correlated or uncorrelated depending on LOCB
    select case(LOCB)
      case(LOCB_CORRELATED)
        self % correlated = .true.
        call new_correlatedLawENDF(self % corrLaw, ACE, MT)

      case default
        self % correlated = .false.
        call new_angleLawENDF(self % muLaw, ACE, MT)
        call new_energyLawENDF(self % eLaw, ACE, MT)

    end select

  end subroutine buildFromACE

  !!
  !! Cast reactionHandle pointer to neutronScatter pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class reactionHandle
  !!
  !! Result:
  !!   Null is source is not of neutronScatter type
  !!   Target points to source if source is neutronScatter type
  !!
  pure function neutronScatter_TptrCast(source) result(ptr)
    class(reactionHandle), pointer, intent(in) :: source
    type(neutronScatter), pointer       :: ptr

    select type(source)
      type is(neutronScatter)
        ptr => source

      class default
        ptr => null()

    end select

  end function neutronScatter_TptrCast

end module neutronScatter_class
