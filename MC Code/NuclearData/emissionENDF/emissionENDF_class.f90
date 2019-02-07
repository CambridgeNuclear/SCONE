module emissionENDF_class

  use numPrecision
  use endfConstants
  use RNG_class,                     only : RNG
  use aceCard_class,                 only : aceCard

  ! Interfaces
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
  !! Constructor
  !!
  interface emissionENDF
    module procedure new_emissionENDF_uncorrelated
    module procedure new_emissionENDF_correlated
    module procedure new_emissionENDF_fromACE
  end interface

  !!
  !! Strores outgoing energy & angle distribution for an MT reaction
  !!
  type, public :: emissionENDF
    private
    ! Reference frame and emission type flags
    logical(defBool)   :: cmFrame    = .true.
    logical(defBool)   :: correlated = .false.

    ! Emission laws
    class(angleLawENDF),allocatable      :: muLaw
    class(energyLawENDF),allocatable     :: eLaw
    class(correlatedLawENDF),allocatable :: corrLaw

  contains
    generic   :: init => init_uncorrelated, init_correlated, init_fromACE

    procedure :: sampleAngleEnergy
    procedure :: isInCMframe
    procedure :: moveAllocFrom

    procedure, private :: init_uncorrelated
    procedure, private :: init_correlated
    procedure, private :: init_fromACE

  end type emissionENDF

contains

  !!
  !! Samples angle mu and energy E_out of outgoing neutron given incident energy and
  !! random number generator
  !!
  subroutine sampleAngleEnergy(self,mu,E_out,E_in,rand )
    class(emissionENDF), intent(in)   :: self
    real(defReal), intent(out)        :: mu
    real(defReal), intent(out)        :: E_out
    real(defReal), intent(in)         :: E_in
    class(RNG), intent(inout)         :: rand

    if(self % correlated) then
      call self % corrLaw % sample(mu,E_out,E_in,rand)

    else
      mu    = self % muLaw % sample(E_in,rand)
      E_out = self % eLaw  % sample(E_in,rand)

    end if

  end subroutine sampleAngleEnergy

  !!
  !! Returns true is data for the reaction is given in Centre-of-Mass Frame
  !!
  function isInCMframe(self)
    class(emissionENDF), intent(in) :: self
    logical(defBool)                :: isInCMframe

    isInCMframe = self % cmFrame
  end function isInCMframe

  !!
  !! Copies RHS into LHS and deallocates contents of RHS
  !!
  subroutine moveAllocFrom(LHS,RHS)
    class(emissionENDF), intent(out)  :: LHS
    type(emissionENDF), intent(inout) :: RHS

    LHS % cmFrame    = RHS % cmFrame
    LHS % correlated = RHS % correlated

    ! Move allocateion to avoid unnecessary memory allocation
    if(RHS % correlated) then
      call move_alloc(RHS % corrLaw, LHS % corrLaw)

    else
      call move_alloc(RHS % muLaw, LHS % muLaw)
      call move_alloc(RHS % eLaw,  LHS % eLaw )

    end if

  end subroutine moveAllocFrom

  !!
  !! Initialisation of uncorrelated emission
  !!
  subroutine init_uncorrelated(self,muLaw,eLaw,cmFrame)
    class(emissionENDF), intent(inout) :: self
    class(angleLawENDF), intent(in)    :: muLaw
    class(energyLawENDF), intent(in)   :: eLaw
    logical(defBool),intent(in)        :: cmFrame

    allocate(self % muLaw, source = muLaw)
    allocate(self % eLaw,  source = eLaw )
    self % correlated = .false.
    self % cmFrame = cmFrame

  end subroutine init_uncorrelated

  !!
  !! Initialisation of correlated emission
  !!
  subroutine init_correlated(self,corrLaw,cmFrame)
    class(emissionENDF), intent(inout)   :: self
    class(correlatedLawENDF), intent(in) :: corrLaw
    logical(defBool),intent(in)          :: cmFrame

    allocate(self % corrLaw, source = corrLaw)
    self % correlated = .true.
    self % cmFrame = cmFrame

  end subroutine init_correlated

  !!
  !! Initialisation of emissionENDF from ACE and MT number
  !! aceCard read head can be in any position
  !! read head will be moved in the subroutine
  !!
  subroutine init_fromACE(self,ACE,MT)
    class(emissionENDF), intent(inout)   :: self
    class(aceCard), intent(inout) :: ACE
    integer(shortInt), intent(in) :: MT
    integer(shortInt)             :: LOCB

    if( ACE % isCaptureMT(MT)) then
      ! Capture Does not have LOCB. Thus build it here.
      ! Will be filled with placeholder LawENDF's
      self % correlated = .false.
      call new_angleLawENDF(self % muLaw, ACE, MT)
      call new_energyLawENDF(self % eLaw, ACE, MT)
      return

    end if

    ! Read LOCB for reaction under MT
    LOCB = ACE % LOCBforMT(MT)

    ! Read if data for reaction is in Centre-of-Mass frame
    self % cmFrame = ACE % isCMframe(MT)

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

  end subroutine init_fromACE

  !!
  !! Constructor of uncorrelated emission
  !!
  function new_emissionENDF_uncorrelated(muLaw,eLaw,cmFrame) result(new)
    class(angleLawENDF), intent(in)  :: muLaw
    class(energyLawENDF), intent(in) :: eLaw
    logical(defBool),intent(in)      :: cmFrame
    type(emissionENDF)               :: new

    call new % init(muLaw,eLaw,cmFrame)

  end function new_emissionENDF_uncorrelated

  !!
  !! Constructor of correlated emission
  !!
  function new_emissionENDF_correlated(corrLaw,cmFrame) result(new)
    class(correlatedLawENDF), intent(in)  :: corrLaw
    logical(defBool),intent(in)           :: cmFrame
    type(emissionENDF)                    :: new

    call new % init(corrLaw,cmFrame)

  end function new_emissionENDF_correlated

  !!
  !! Constructor from ACE and MT number
  !! aceCard read head can be in any position
  !! read head will be moved in the subroutine
  !!
  function new_emissionENDF_fromACE(ACE,MT) result(new)
    class(aceCard), intent(inout) :: ACE
    integer(shortInt), intent(in) :: MT
    type(emissionENDF)            :: new

    call new % init(ACE,MT)

  end function new_emissionENDF_fromACE

end module emissionENDF_class
