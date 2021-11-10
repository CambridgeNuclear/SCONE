module thermalScatterElastic_class

  use numPrecision
  use endfConstants
  use genericProcedures,            only : fatalError, numToChar, binarySearch
  use RNG_class,                    only : RNG
  use dataDeck_inter,               only : dataDeck
  use aceSabCard_class,             only : aceSabCard
  use reactionHandle_inter,         only : reactionHandle
  use uncorrelatedReactionCE_inter, only : uncorrelatedReactionCE

  implicit none
  private

  !!
  !! Public pointer cast
  !!
  public :: thermalScatterElastic_TptrCast

  !!
  !!
  !!
  type, public, extends(uncorrelatedReactionCE) :: thElasticScatter
    private
    real(defReal), dimension(:), allocatable   :: eIn
    real(defReal), dimension(:,:), allocatable :: muMatrix
    logical(defBool)   :: elasticTable
    integer(shortInt)  :: N_muOut
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

  end type thElasticScatter

contains

  !!
  !! Initialise
  !!
  !! See reactionHandle for details
  !!
  !! Errors:
  !!   fatalError if MT is not neutron elastic scattering
  !!
  subroutine init(self, data, MT)
    class(thElasticScatter), intent(inout) :: self
    class(dataDeck), intent(inout)         :: data
    integer(shortInt), intent(in)          :: MT
    character(100), parameter :: Here = 'init (elasticNeutronScatter_class.f90)'

    ! Select buld procedure approperiate for given dataDeck
    select type(data)
      type is (aceSabCard)
        call self % buildFromACE(data)

      class default
        call fatalError(Here,'Elastic thermal neutron scattering cannot be build from '//data % myType())
    end select

  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  !! See reactionHandle for details
  !!
  elemental subroutine kill(self)
    class(thElasticScatter), intent(inout) :: self

  end subroutine kill

  !!
  !! Returns true if reaction is in Centre-Of-Mass frame
  !!
  !! See uncorrelatedReactionCE for details
  !!
  pure function inCMframe(self) result(isIt)
    class(thElasticScatter), intent(in) :: self
    logical(defBool)                    :: isIt

     isIt = .true.

  end function inCMframe

  !!
  !! Returns number of particles produced on average by the reaction
  !!
  !! See uncorrelatedReactionCE for details
  !!
  pure function release(self, E) result(N)
    class(thElasticScatter), intent(in) :: self
    real(defReal), intent(in)           :: E
    real(defReal)                       :: N

    N = ONE

  end function release

  !!
  !! Returns number of particles produced on average instantly by the reaction
  !!
  !! See uncorrelatedReactionCE for details
  !!
  pure function releasePrompt(self, E) result(N)
    class(thElasticScatter), intent(in) :: self
    real(defReal), intent(in)           :: E
    real(defReal)                       :: N

    N = ONE

  end function releasePrompt

  !!
  !! Returns number of particles produced on average by the reaction with delay
  !!
  !! See uncorrelatedReactionCE for details
  !!
  pure function releaseDelayed(self, E) result(N)
    class(thElasticScatter), intent(in) :: self
    real(defReal), intent(in)           :: E
    real(defReal)                       :: N

    N = ZERO

  end function releaseDelayed

  !!
  !! Sample outgoing particle
  !!
  !! See uncorrelatedReactionCE for details
  !!
  subroutine sampleOut(self, mu, phi, E_out, E_in, rand, lambda)
    class(thElasticScatter), intent(in) :: self
    real(defReal), intent(out)               :: mu
    real(defReal), intent(out)               :: phi
    real(defReal), intent(out)               :: E_out
    real(defReal), intent(in)                :: E_in
    class(RNG), intent(inout)                :: rand
    real(defReal), intent(out), optional     :: lambda
    real(defReal)        :: E1, E2, rho, mu_lk, mu1, mu2, mu3, muLeft, muRight
    integer(shortInt)    :: l1, l2, k
    character(100), parameter :: Here = 'sampleOut(thermalScatterElastic_class)'

    ! Set energy
    E_out = E_in

    if (self % elasticTable .and. E_in < self % eIn(size(self % eIn))) then

      ! Get energy indexes
      l1 = binarySearch(self % eIn, E_in)
      l2 = l1 + 1

      ! Get energy values
      E1 = self % eIn(l1)
      E2 = self % eIn(l2)
      rho = (E_in - E1)/(E2 - E1)

      ! Sample the outgoing angle
      k = floor(self % N_muOut * rand % get()) + 1
      mu_lk = self % muMatrix(l1,k)
      mu1 = mu_lk + rho * (self % muMatrix(l2,k) - mu_lk)

      ! Smearing the outgoing angular distribution
      if (k /= 1 .and. k /= self % N_muOut) then
        mu2 = self % muMatrix(l1, k - 1)
        mu3 = self % muMatrix(l2, k - 1)
        muLeft = mu2 + rho * (mu3 - mu2)
        mu2 = self % muMatrix(l1, k + 1)
        mu3 = self % muMatrix(l2, k + 1)
        muRight = mu2 + rho * (mu3 - mu2)
        mu = mu1 + min((mu1 - muLeft),(mu1 + muRight))*(rand % get() - 0.5_defReal)
        if (mu > ONE .or. mu < -ONE) mu = mu1
      else
        mu = mu1
      end if

    else
      if (E_in >= self % eIn(size(self % eIn))) then
        E1 = self % eIn(size(self % eIn))
      else
        ! Get energy indexes
        l1 = binarySearch(self % eIn, E_in)
        E1 = self % eIn(l1)
      end if
      mu = 1 - E1/E_in

    end if

    if (mu > ONE .or. mu < -ONE) then
      call fatalError(Here,'Failed to get angle'//numToChar(mu))
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
    class(thElasticScatter), intent(in) :: self
    real(defReal), intent(in)           :: mu
    real(defReal), intent(in)           :: phi
    real(defReal), intent(in)           :: E_out
    real(defReal), intent(in)           :: E_in
    real(defReal)                       :: prob

    prob = ONE

  end function probOf

  !!
  !! Build elasticNeutronScatter from ACE dataCard
  !!
  subroutine buildFromACE(self, ACE)
    class(thElasticScatter), intent(inout) :: self
    type(aceSabCard), intent(inout)        :: ACE
    integer(shortInt)                      :: Nin, i
    character(100), parameter :: Here = 'buildFromACE (thElasticScatt_class.f90)'

    Nin = ACE % inelasticEnergies()
    self % eIn = ACE % ESZ_elastic('energyGrid')

    self % elasticTable = ACE % hasITCA()

    if (self % elasticTable) then

      self % N_muOut = ACE % elOutMu()
      call ACE % setToElasticOut()
      allocate(self % muMatrix(Nin, self % N_muOut))

      do i = 1, Nin
        self % muMatrix(i,:) = ACE % readIntArray(Nin)
      end do

    end if

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
  pure function thermalScatterElastic_TptrCast(source) result(ptr)
    class(reactionHandle), pointer, intent(in) :: source
    type(thElasticScatter), pointer            :: ptr

    select type(source)
      type is(thElasticScatter)
        ptr => source

      class default
        ptr => null()
    end select

  end function thermalScatterElastic_TptrCast


end module thermalScatterElastic_class
