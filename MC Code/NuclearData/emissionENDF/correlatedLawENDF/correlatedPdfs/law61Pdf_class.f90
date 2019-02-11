module law61Pdf_class

  use numPrecision
  use endfConstants
  use genericProcedures,   only : fatalError
  use RNG_class,           only : RNG
  use aceCard_class,       only : aceCard
  use muEndfPdfSlot_class, only : muEndfPdfSlot
  use tabularEnergy_class, only : tabularEnergy

  implicit none
  private

  !!
  !! This pdf stores an outgoing energy grid and outgoing angle-pdfs associated with each point
  !! Each pdf corresponds to a single point on incident energy grid in ACE LAW 61 formulation
  !!
  type, public :: law61Pdf
    private
    type(tabularEnergy)                           :: ePdf
    type(muEndfPdfSlot),dimension(:), allocatable :: muPdfs
  contains
    !! Public Interface
    generic   :: init => init_fromACE
    procedure :: sample
    procedure :: probabilityOf
    procedure :: bounds

    !! Private procedures
    procedure, private :: init_fromACE
  end type law61Pdf

contains

  !!
  !! Samples mu and E_out given random nummber generator
  !!
  subroutine sample(self, mu, E_out, rand)
    class(law61Pdf), intent(in) :: self
    real(defReal), intent(out)  :: mu
    real(defReal), intent(out)  :: E_out
    class(RNG), intent(inout)   :: rand
    real(defReal)               :: eps, r
    integer(shortInt)           :: bin

    ! Sample Energy
    E_out = self % ePdf % sample(rand, bin)

    ! Sample Angle
    eps = self % ePdf % getInterF(E_out, bin)
    r = rand % get()

    if(r < eps) then
      mu = self % muPdfs(bin+1) % sample(rand)
    else
      mu = self % muPdfs(bin) % sample(rand)

    end if

  end subroutine sample

  !!
  !! Returns probability that neutron was emmited at mu & E_out
  !!
  function probabilityOf(self, mu, E_out) result(prob)
    class(law61Pdf), intent(in) :: self
    real(defReal), intent(in)   :: mu
    real(defReal), intent(in)   :: E_out
    real(defReal)               :: prob
    real(defReal)               :: f
    integer(shortInt)           :: bin

    ! Probability of a given energy
    prob = self % ePdf % probabilityOf(E_out, bin)

    ! Propability of a given angle at the energy
    f = self % ePdf % getInterF(E_out, bin)
    prob = prob * ( (ONE-f) * self % muPdfs(bin) % probabilityOf(mu) + &
                          f * self % muPdfs(bin+1) % probabilityOf(mu))

  end function probabilityOf

  !!
  !! Return energy bounds of the probability distribution
  !!
  subroutine bounds(self, E_min, E_max)
    class(law61Pdf), intent(in) :: self
    real(defReal), intent(out)  :: E_min
    real(defReal), intent(out)  :: E_max

    call self % ePdf % bounds(E_min, E_max)

  end subroutine bounds

  !!
  !! Initialise from ACE. Overwrite and extend superclass procedure
  !!
  subroutine init_fromACE(self, ACE)
    class(law61Pdf), intent(inout)              :: self
    class(aceCard), intent(inout)               :: ACE
    integer(shortInt)                           :: NP, i
    integer(shortInt),dimension(:), allocatable :: LCs
    character(100), parameter :: Here = 'init_fromACE (law61Pdf_class.f90)'

    ! Read number of points without moving read head
    call ACE % advanceHead(1)
    NP = ACE % readIntNotAdvance()
    call ACE % advanceHEad(-1)

    ! Call superclass initialisation
    call self % ePdf % init(ACE)

    ! Read LC table (locators of angular distributions) and allocate space
    LCs = ACE % readIntArray(NP)
    allocate(self % muPdfs(NP))

    ! Read angular data
    do i=1,NP
      select case(LCs(i))
        case(1:)
          call ACE % setToEnergyLaw(LCs(i))
          call self % muPdfs(i) % init(ACE, 'tabularMu')

        case (0)
          call self % muPdfs(i) % init(ACE, 'isotropicMu')

        case(:-1)
          call fatalError(Here,'32 Equiprobable bin angular distribution to allowed in LAW 61!')

      end select
    end do

  end subroutine init_fromACE



end module law61Pdf_class
