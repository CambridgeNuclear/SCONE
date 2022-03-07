module thermalScatterInelastic_class

  use numPrecision
  use endfConstants
  use genericProcedures,            only : binarySearch, fatalError, numToChar, endfInterpolate, searchError
  use RNG_class,                    only : RNG
  use dataDeck_inter,               only : dataDeck
  use aceSabCard_class,             only : aceSabCard
  use tabularEnergy_class,          only : tabularEnergy
  use reactionHandle_inter,         only : reactionHandle
  use uncorrelatedReactionCE_inter, only : uncorrelatedReactionCE

  implicit none
  private

  !!
  !! Public pointer cast
  !!
  public :: thermalScatterInelastic_TptrCast

  !! Private type to store angular tables
  type, private :: angularMatrix
    real(defReal), dimension(:,:), allocatable :: muOut
  end type angularMatrix

  !! Private type to store energy and CDF vectors
  type, private :: dataArray
    real(defReal), dimension(:), allocatable :: array
  end type dataArray

  !!
  !! Reaction type for Neutron Elastic Scattering
  !!
  !! Implements standard elastic neutron stattering
  !! For conveniance can exist in two states: isotropic & anisotropic.
  !! Data is assumed to be always in CoM-frame.
  !!
  !! NOTE:
  !!   When building from ACE data don't be fooled by Appendix F of MCNP-4 manual.
  !!   It is possible for LOCB for elastic scattering to be 0, which indicates isotropic
  !!   scattering at all incident energies. Thus LOCB needs to be checked and `isotropic` flag
  !!   set to a correct setting.
  !!
  !! Private Members:
  !!   angularData -> tabularAngle type that holds the energy-dependant data
  !!   isotropic   -> Flag, TRUE is scattering is isotropic at all incident energies
  !!
  !! Interface:
  !!   uncorrelatedReactionCE interface
  !!   buildFromACE -> initialise object from ACE dataCard
  !!
  type, public, extends(uncorrelatedReactionCE) :: thInelasticScatter
    private
    real(defReal), dimension(:), allocatable       :: eIn
    real(defReal), dimension(:), allocatable       :: CDF
    type(tabularEnergy),dimension(:),allocatable   :: eOutPdf
    type(dataArray), dimension(:), allocatable     :: eOut
    type(angularMatrix), dimension(:), allocatable :: muMatrices
    logical(defBool)   :: isInelContinuous
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

  end type thInelasticScatter

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
    class(thInelasticScatter), intent(inout) :: self
    class(dataDeck), intent(inout)           :: data
    integer(shortInt), intent(in)            :: MT
    character(100), parameter :: Here = 'init (thermalScatterInelastic_class.f90)'

    ! Select build procedure appropriate for given dataDeck
    select type(data)
    type is (aceSabCard)
        call self % buildFromACE(data)
      class default
        call fatalError(Here,'Inelastic thermal neutron scattering cannot be build from '//data % myType())
    end select

  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  !! See reactionHandle for details
  !!
  elemental subroutine kill(self)
    class(thInelasticScatter), intent(inout) :: self

  end subroutine kill

  !!
  !! Returns true if reaction is in Centre-Of-Mass frame
  !!
  !! See uncorrelatedReactionCE for details
  !!
  pure function inCMframe(self) result(isIt)
    class(thInelasticScatter), intent(in) :: self
    logical(defBool)                         :: isIt

     isIt = .false.

  end function inCMframe

  !!
  !! Returns number of particles produced on average by the reaction
  !!
  !! See uncorrelatedReactionCE for details
  !!
  pure function release(self, E) result(N)
    class(thInelasticScatter), intent(in) :: self
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
    class(thInelasticScatter), intent(in) :: self
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
    class(thInelasticScatter), intent(in) :: self
    real(defReal), intent(in)           :: E
    real(defReal)                       :: N

    N = ZERO

  end function releaseDelayed

  !!
  !! Return probability density of emission at given angle and energy
  !!
  !! See uncorrelatedReactionCE for details
  !!
  function probOf(self, mu, phi, E_out, E_in) result(prob)
    class(thInelasticScatter), intent(in) :: self
    real(defReal), intent(in)             :: mu
    real(defReal), intent(in)             :: phi
    real(defReal), intent(in)             :: E_out
    real(defReal), intent(in)             :: E_in
    real(defReal)                         :: prob

    prob = ZERO


  end function probOf

  !!
  !! Sample outgoing particle
  !!
  !! See uncorrelatedReactionCE for details
  !!
  subroutine sampleOut(self, mu, phi, E_out, E_in, rand, lambda)
    class(thInelasticScatter), intent(in) :: self
    real(defReal), intent(out)            :: mu
    real(defReal), intent(out)            :: phi
    real(defReal), intent(out)            :: E_out
    real(defReal), intent(in)             :: E_in
    class(RNG), intent(inout)             :: rand
    real(defReal), intent(out), optional  :: lambda
    real(defReal)     :: E_min_low, E_max_low
    real(defReal)     :: E_min_up, E_max_up
    real(defReal)     :: E_min, E_max
    real(defReal)     :: E1, E2, f, ri, eps
    real(defReal)     :: mu_ljk, mu1, mu2, mu3, muLeft, muRight
    integer(shortInt) :: l1, l2, l, j, k
    character(100), parameter :: Here = 'sampleOut(thermalScatterInelastic_class)'

    ! Get energy indexes
    l1 = binarySearch(self % eIn, E_in)
    l2 = l1 + 1

    ! Get energy values
    E1 = self % eIn(l1)
    E2 = self % eIn(l2)

    f = (E_in - E1)/(E2 - E1)

    if ( .not. self % isInelContinuous) then
      j = binarySearch(self % CDF, rand % get())
      E_min = self % eOut(l1) % array(j)
      E_max = self % eOut(l2) % array(j)

      E_out = E_min + f * (E_max - E_min)

      ! Considering a uniform distribution for the angular bins
      k = floor(self % N_muOut * rand % get()) + 1
      mu2 = self % muMatrices(l1) % muOut(j, k)
      mu3 = self % muMatrices(l2) % muOut(j, k)

      mu = mu2 + f * (mu3 - mu2)

    else
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculate energy boundaries
      call self % eOutPdf(l1) % bounds(E_min_low, E_max_low)
      call self % eOutPdf(l2) % bounds(E_min_up,  E_max_up )

      ! Calculate interpolated bounds
      E_min = E_min_low * (ONE - f) + f * E_min_up
      E_max = E_max_low * (ONE - f) + f * E_max_up

      if(rand % get() < f) then
        l = l2
        E_out = self % eOutPdf(l) % sample(rand, j)
        ri = (E_out - E_min_up)/(E_max_up - E_min_up)
      else
        l = l1
        E_out = self % eOutPdf(l) % sample(rand, j)
        ri = (E_out - E_min_low)/(E_max_low - E_min_low)
      end if

      eps = self % eOutPdf(l) % getInterF(E_out, j)
      E_out = E_min + ri*(E_max - E_min)

      ! Sampling the outgoing angle
      k = floor(self % N_muOut * rand % get()) + 1
      mu_ljk = self % muMatrices(l) % muOut(j, k)
      mu1 = mu_ljk + eps * (self % muMatrices(l) % muOut(j + 1, k) - mu_ljk)

      ! Smearing the outgoing angular distribution
!      if (k /= 1 .and. k /= self % N_muOut) then
!        mu2 = self % muMatrices(l) % muOut(j, k - 1)
!        mu3 = self % muMatrices(l) % muOut(j + 1, k - 1)
!        muLeft = mu2 + eps * (mu3 - mu2)
!        mu2 = self % muMatrices(l) % muOut(j, k + 1)
!        mu3 = self % muMatrices(l) % muOut(j + 3, k + 1)
!        muRight = mu2 + eps * (mu3 - mu2)
!        mu = mu1 + min((mu1 - muLeft),(mu1 + muRight))*(rand % get() - 0.5_defReal)
!        if (mu > ONE .or. mu < -ONE) mu = mu1
!      else
!        mu = mu1
!      end if
      mu = mu1
      if (mu > ONE .or. mu < -ONE) mu = mu_ljk
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end if

    if (E_out > 20.0 .or. E_out <= ZERO) then
      call fatalError(Here,'Failed to find energy: '//numToChar(E_out))
    end if

    if (.not. E_out == E_out ) print*, E_out, E_min, E_max, ri, j

    if (mu > ONE .or. mu < -ONE) then
      print*, self % N_muOut, k, j, eps, mu_ljk, mu1, mu2, mu3, muLeft, muRight
      call fatalError(Here,'Failed to get angle'//numToChar(mu))
    end if

    ! Sample phi
    phi = rand % get() * TWO_PI

    ! Only prompt particles. Set delay
    if(present(lambda)) lambda = huge(lambda)

  end subroutine sampleOut

  !!
  !! Build elasticNeutronScatter from ACE dataCard
  !!
  subroutine buildFromACE(self, ACE)
    class(thInelasticScatter), intent(inout)      :: self
    type(aceSabCard), intent(inout)               :: ACE
    real(defReal), dimension(:), allocatable      :: Etmp, PDFtmp, CDFtmp
    integer(shortInt), dimension(:), allocatable  :: loc, Eout
    integer(shortInt)                             :: i, j, Nin, Nout
    character(100), parameter :: Here = 'buildFromACE (thInelasticScatt_class.f90)'

    Nin = ACE % inelasticEnergies()
    self % eIn = ACE % ESZ_inelastic('energyGrid')
    self % isInelContinuous = ACE % isInelContinuous()

    call ACE % setToInelasticOut()

    if (.not. self % isInelContinuous) then
      self % N_muOut = ACE % inelOutMu() + 1
      Nout = ACE % inelOutE()

      allocate(self % eOut(Nin), self % muMatrices(Nin), &
                Etmp(Nout), self % CDF(Nout + 1))
      self % CDF = ACE % getCDF()

      do i = 1, Nin
        allocate(self % muMatrices(i) % muOut(Nout, self % N_muOut))
        do j = 1, Nout
          Etmp(j) = ACE % readReal()
          self % muMatrices(i) % muOut(j,:) = ACE % readRealArray(self % N_muOut)
        end do
        allocate(self % eOut(i) % array(Nout))
        self % eOut(i) % array = Etmp(1:Nout)
      end do

    else

      self % N_muOut = ACE % inelOutMu() - 1
      allocate(self % eOutPdf(Nin), self % muMatrices(Nin))

      loc = ACE % readIntArray(Nin)
      Eout = ACE % readIntArray(Nin)

      do i = 1, Nin

        allocate(self % muMatrices(i) % muOut(Eout(i), self % N_muOut), &
        Etmp(Eout(i)), PDFtmp(Eout(i)), CDFtmp(Eout(i)+1))

        call ACE % setRelativeTo(loc(i), 1)
        CDFtmp(1) = ZERO
        do j = 1, Eout(i)
          Etmp(j) = ACE % readReal()
          PDFtmp(j) = ACE % readReal()
          CDFtmp(j+1) = ACE % readReal()
          self % muMatrices(i) % muOut(j,:) = ACE % readRealArray(self % N_muOut)
        end do

        PDFtmp(Eout(i)) = ZERO
        self % eOutPdf(i) = tabularEnergy(Etmp, PDFtmp, CDFtmp(1:Eout(i)), 2)

        deallocate(Etmp, PDFtmp, CDFtmp)

      end do

    end if

  end subroutine buildFromACE

  !!
  !! Cast reactionHandle pointer to thInelasticScatt pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class reactionHandle
  !!
  !! Result:
  !!   Null is source is not of elasticNeutronScatter type
  !!   Target points to source if source is elasticNeutronScatter type
  !!
  pure function thermalScatterInelastic_TptrCast(source) result(ptr)
    class(reactionHandle), pointer, intent(in) :: source
    type(thInelasticScatter), pointer          :: ptr

    select type(source)
      type is(thInelasticScatter)
        ptr => source

      class default
        ptr => null()
    end select

  end function thermalScatterInelastic_TptrCast


end module thermalScatterInelastic_class
