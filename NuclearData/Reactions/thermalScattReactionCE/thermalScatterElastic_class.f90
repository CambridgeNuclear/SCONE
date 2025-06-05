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
  !! Reaction type for Neutron Thermal Elastic Scattering with S(a,b) tables
  !!
  !! The outgoing energy is always the same as the ingong energy. The angles can
  !! be either tabulated or calculated after sampling a Bragg edge:
  !! https://docs.openmc.org/en/stable/methods/neutron_physics.html#outgoing-angle-for-coherent-elastic-scattering
  !!
  !! Data is always in lab-frame.
  !!
  !! Private Members:
  !!   eIn          -> ingoing energy grid
  !!   pValues      -> values from ACE file - could be xss or data to be processed to obtain xss
  !!   muMatrix     -> outgoing angles matrix
  !!   elasticTable -> flag about whether angular data are tabulate or not
  !!   N_muOut      -> number of outgoing angle cosines
  !!
  !! Interface:
  !!   uncorrelatedReactionCE interface
  !!   buildFromACE -> initialise object from ACE Sab dataCard
  !!
  type, public, extends(uncorrelatedReactionCE) :: thElasticScatter
    private
    real(defReal), dimension(:), allocatable   :: eIn
    real(defReal), dimension(:), allocatable   :: pValues
    real(defReal), dimension(:,:), allocatable :: muMatrix
    logical(defBool)   :: elasticTable = .false.
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

    self % N_muOut = 0
    self % elasticTable = .false.

    if(allocated(self % muMatrix)) deallocate(self % muMatrix)
    deallocate(self % eIn)
    deallocate(self % pValues)

  end subroutine kill

  !!
  !! Returns true if reaction is in Centre-Of-Mass frame
  !!
  !! See uncorrelatedReactionCE for details
  !!
  pure function inCMframe(self) result(isIt)
    class(thElasticScatter), intent(in) :: self
    logical(defBool)                    :: isIt

     isIt = .false.

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
    real(defReal), dimension(:), allocatable :: prob
    real(defReal)        :: E1, E2, f, mu_l1k, mu1, mu2, mu3, muLeft, muRight, r
    integer(shortInt)    :: l1, l2, k, i
    character(100), parameter :: Here = 'sampleOut(thermalScatterElastic_class)'

    ! Set energy
    E_out = E_in

    if (self % elasticTable) then

      ! Get energy indexes
      if (E_in > self % eIn(size(self % eIn))) then
        l1 = size(self % eIn) - 1
      else
        l1 = binarySearch(self % eIn, E_in)
      end if

      l2 = l1 + 1

      ! Get energy values
      E1 = self % eIn(l1)
      E2 = self % eIn(l2)
      f = (E_in - E1)/(E2 - E1)

      ! Sampling loop
      sample: do i = 1, 100

        ! Sample the outgoing angle
        k = floor(self % N_muOut * rand % get()) + 1
        mu_l1k = self % muMatrix(l1,k)
        mu1 = mu_l1k + f * (self % muMatrix(l2,k) - mu_l1k)

        ! Smearing the outgoing angular distribution
        if (k == 1) then
          muLeft = - ONE - (mu1 + ONE)
        else
          mu2 = self % muMatrix(l1, k - 1)
          mu3 = self % muMatrix(l2, k - 1)
          muLeft = mu2 + f * (mu3 - mu2)
        end if

        if (k == self % N_muOut) then
          muRight = ONE + (ONE - mu1)
        else
          mu2 = self % muMatrix(l1, k + 1)
          mu3 = self % muMatrix(l2, k + 1)
          muRight = mu2 + f * (mu3 - mu2)
        end if

        mu = mu1 + min(mu1 - muLeft, muRight - mu1) * (rand % get() - HALF)

        ! Check if the angle is valid
        if (mu <= ONE .and. mu >= - ONE) exit sample

      end do sample

      if (i == 100) call fatalError(Here,'Failed to find angle: '//numToChar(mu))

    else

      ! Get energy indexes
      if (E_in > self % eIn(size(self % eIn))) then
        l1 = size(self % eIn) - 1
      else
        l1 = binarySearch(self % eIn, E_in)
      end if

      l2 = l1 + 1

      ! Sample a Bragg edge at a smaller energy than E_in
      prob = self % pValues(1:l2)
      prob(1) = ZERO
      r = rand % get() * prob(l2)
      k = binarySearch(prob, r)

      E2 = self % eIn(k)
      ! Compute angle
      mu = 1 - TWO * E2/E_in

      if (abs(mu) > ONE) then
        call fatalError(Here,'Failed to get angle'//numToChar(mu))
      end if
    end if

    ! Sample phi
    phi = rand % get() * TWO_PI

    ! Only prompt particles. Set delay
    if(present(lambda)) lambda = huge(lambda)

  end subroutine sampleOut

  !!
  !! Return probability density of emission at given angle and energy
  !! NOT IMPLEMENTED: raises fatal error if called
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
    character(100), parameter :: Here = 'probOf (thermalScatterElastic_class.f90)'

    ! Avoid compiler warnings
    prob = ONE
    ! Not implemented yet
    call fatalError(Here,'This function is not implemented')

  end function probOf

  !!
  !! Build elasticNeutronScatter from ACE dataCard
  !!
  subroutine buildFromACE(self, ACE)
    class(thElasticScatter), intent(inout) :: self
    type(aceSabCard), intent(inout)        :: ACE
    integer(shortInt)                      :: Nin, i
    character(100), parameter :: Here = 'buildFromACE (thElasticScatt_class.f90)'

    ! Get scattering grids
    self % eIn = ACE % ESZ_elastic('energyGrid')
    self % Pvalues = ACE % ESZ_elastic('Pvalues')

    ! Check if outgoing distribution is tabulated
    self % elasticTable = ACE % hasITCA()

    ! Get tabulated outgoing angle distributions
    if (self % elasticTable) then

      Nin = ACE % elasticEnergies()
      self % N_muOut = ACE % elOutMu()
      allocate(self % muMatrix(Nin, self % N_muOut))

      call ACE % setToElasticOut()
      do i = 1, Nin
        self % muMatrix(i,:) = ACE % readRealArray(self % N_muOut)
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
