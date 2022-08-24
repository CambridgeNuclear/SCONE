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
  !! Reaction type for Neutron Thermal Inelastic Scattering with S(a,b) tables
  !!
  !! Implements two outgoing scattering laws: discrete and continuous tables. Continuous
  !! tables are implemented from library ENDF/B-VIII.0 on.
  !!
  !! Data is always in lab-frame.
  !!
  !! Private Members:
  !!   eIn        -> ingoing energy grid
  !!   CDF        -> CDF for outgoing energy distribution
  !!   eOutPdf    -> continuous tabular energy endf law
  !!   eOut       -> array of outgoing energies arrays (which might have different lengths)
  !!   muMatrices -> array of outgoing angles/energies matrices (which might have different sizes)
  !!   isInelContinuous -> flag about continuous tabular format
  !!   N_muOut    -> number of cosine bins
  !!
  !! Interface:
  !!   uncorrelatedReactionCE interface
  !!   buildFromACE -> initialise object from ACE Sab dataCard
  !!
  type, public, extends(uncorrelatedReactionCE) :: thInelasticScatter
    private
    real(defReal), dimension(:), allocatable       :: eIn
    real(defReal), dimension(:), allocatable       :: CDF
    type(tabularEnergy),dimension(:),allocatable   :: eOutPdf
    type(dataArray), dimension(:), allocatable     :: eOut
    type(angularMatrix), dimension(:), allocatable :: muMatrices
    logical(defBool)   :: isInelContinuous = .false.
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

    self % N_muOut = 0
    self % isInelContinuous = .false.

    if(allocated(self % eOut))    deallocate(self % eOut)
    if(allocated(self % CDF))     deallocate(self % CDF)
    if(allocated(self % eOutPdf)) deallocate(self % eOutPdf)
    deallocate(self % eIn)
    deallocate(self % muMatrices)

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
  !! NOT IMPLEMENTED: raises fatal error if called
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
    character(100), parameter :: Here = 'probOf (thermalScatterInelastic_class.f90)'

    ! Avoid compiler warnings
    prob = ONE
    ! Not implemented yet
    call fatalError(Here,'This function is not implemented')

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
    real(defReal)     :: E_min, E_max
    real(defReal)     :: E1, E2, f, eps
    real(defReal)     :: mu_ljk, mu1, mu2, mu3, muLeft, muRight
    integer(shortInt) :: l1, l2, l, j, k, i
    character(100), parameter :: Here = 'sampleOut(thermalScatterInelastic_class)'

    ! Get energy indexes
    l1 = binarySearch(self % eIn, E_in)
    l2 = l1 + 1

    ! Get energy values
    E1 = self % eIn(l1)
    E2 = self % eIn(l2)

    f = (E_in - E1)/(E2 - E1)

    if ( .not. self % isInelContinuous) then
      ! Discrete treatment
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
      ! Continuous treatment
      ! Choose index from ingoing energy grid
      if (f > HALF) then
        l = l2
      else
        l = l1
      end if

      ! Sampling loop
      sample: do i = 1, 100

        ! Sample outgoing energy
        E_out = self % eOutPdf(l) % sample(rand, j, eps)

        ! Adjust outgoing energy
        if (E_out < self % eIn(l) * HALF) then
          E_out = E_out * TWO * E_in / self % eIn(l) - E_out
        else
          E_out = E_out + E_in - self % eIn(l)
        end if

        ! Sampling the outgoing angle
        k = floor(self % N_muOut * rand % get()) + 1
        mu_ljk = self % muMatrices(l) % muOut(j, k)
        mu1 = mu_ljk + eps * (self % muMatrices(l) % muOut(j + 1, k) - mu_ljk)

        ! Smearing the outgoing angular distribution
        if (k == 1) then
          muLeft = - ONE - (mu1 + ONE)
        else
          mu2 = self % muMatrices(l) % muOut(j, k - 1)
          mu3 = self % muMatrices(l) % muOut(j + 1, k - 1)
          muLeft = mu2 + eps * (mu3 - mu2)
        end if

        if (k == self % N_muOut) then
          muRight = ONE + (ONE - mu1)
        else
          mu2 = self % muMatrices(l) % muOut(j, k + 1)
          mu3 = self % muMatrices(l) % muOut(j + 1, k + 1)
          muRight = mu2 + eps * (mu3 - mu2)
        end if

        mu = mu1 + min(mu1 - muLeft, muRight - mu1) * (rand % get() - HALF)

        ! Check if the angle is valid
        if (mu <= ONE .and. mu >= - ONE) exit sample

      end do sample

      if (i == 100) call fatalError(Here,'Failed to find angle: '//numToChar(mu))

    end if

    if (E_out > 20.0 .or. E_out <= ZERO) then
      call fatalError(Here,'Failed to find energy: '//numToChar(E_out))
    end if

    ! Sample phi
    phi = rand % get() * TWO_PI

    ! Only prompt particles. Set delay
    if(present(lambda)) lambda = huge(lambda)

  end subroutine sampleOut

  !!
  !! Build inelasticNeutronScatter from ACE dataCard
  !!
  !! See uncorrelatedReactionCE for details
  !!
  subroutine buildFromACE(self, ACE)
    class(thInelasticScatter), intent(inout)      :: self
    type(aceSabCard), intent(inout)               :: ACE
    real(defReal), dimension(:), allocatable      :: Etmp, PDFtmp, CDFtmp
    integer(shortInt), dimension(:), allocatable  :: loc, Eout
    integer(shortInt)                             :: i, j, Nin, Nout
    character(100), parameter :: Here = 'buildFromACE (thInelasticScatt_class.f90)'

    ! Initialise flags from ACE file
    Nin = ACE % inelasticEnergies()
    self % eIn = ACE % ESZ_inelastic('energyGrid')
    self % isInelContinuous = ACE % isInelContinuous()

    call ACE % setToInelasticOut()

    ! Initialise distributions for discrete format
    if (.not. self % isInelContinuous) then
      self % N_muOut = ACE % inelOutMu() + 1
      Nout = ACE % inelOutE()

      allocate(self % eOut(Nin), self % muMatrices(Nin), &
                Etmp(Nout), self % CDF(Nout + 1))

      ! Get CDF distribution built in the ACE file
      self % CDF = ACE % getCDF()
      ! Loop over ingoing energies
      do i = 1, Nin
        allocate(self % muMatrices(i) % muOut(Nout, self % N_muOut))
        ! Loop over outgoing energies
        do j = 1, Nout
          Etmp(j) = ACE % readReal()
          self % muMatrices(i) % muOut(j,:) = ACE % readRealArray(self % N_muOut)
        end do
        allocate(self % eOut(i) % array(Nout))
        self % eOut(i) % array = Etmp(1:Nout)
      end do

    else ! Initialise continuous treatment

      self % N_muOut = ACE % inelOutMu() - 1
      allocate(self % eOutPdf(Nin), self % muMatrices(Nin))

      loc = ACE % readIntArray(Nin)
      Eout = ACE % readIntArray(Nin)
      ! Loop over ingoing energies
      do i = 1, Nin

        allocate(self % muMatrices(i) % muOut(Eout(i) + 1, self % N_muOut), &
        Etmp(Eout(i)+1), PDFtmp(Eout(i)+1), CDFtmp(Eout(i)+1))

        call ACE % setRelativeTo(loc(i), 1)
        ! Make sure initial values are zero
        Etmp(1) = ZERO
        PDFtmp(1) = ZERO
        CDFtmp(1) = ZERO
        self % muMatrices(i) % muOut(1,:) = ZERO
        ! loop over outgoing energies
        do j = 2, Eout(i) + 1
          Etmp(j) = ACE % readReal()
          PDFtmp(j) = ACE % readReal()
          CDFtmp(j) = ACE % readReal()
          self % muMatrices(i) % muOut(j,:) = ACE % readRealArray(self % N_muOut)
        end do

        self % muMatrices(i) % muOut(1,:) = self % muMatrices(i) % muOut(2,:)
        PDFtmp(Eout(i)+1) = ZERO
        ! Use tabular energy class, which includes checks on CDF distribution
        self % eOutPdf(i) = tabularEnergy(Etmp, PDFtmp, CDFtmp, 2)

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
  !!   Null is source is not of inelasticNeutronScatter type
  !!   Target points to source if source is inelasticNeutronScatter type
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
