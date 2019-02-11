module tabularEnergy_class

  use numPrecision
  use genericProcedures,   only : fatalError
  use aceCard_class,       only : aceCard
  use tabularPdf_class,    only : tabularPdf
  use RNG_class,           only : RNG

  implicit none
  private

  interface tabularEnergy
    module procedure new_tabularEnergy
    module procedure new_tabularEnergy_withCDF
    module procedure new_tabularEnergy_fromACE
  end interface

  !!
  !! Probability distribuition table for energy of emitted neutron at a single incomeing energy
  !! Does not support discrete photon lines
  !! NOTE: It is effectivly a wrapper around tabularPdf with extra checks for validity of
  !!       the data (Energies must be +ve)
  !!
  type, public:: tabularEnergy
    private
    type(tabularPdf) :: pdf
  contains
    procedure :: sample
    procedure :: bounds
    procedure :: probabilityOf

    generic   :: init          => init_withPDF, init_withCDF, init_fromACE
    generic   :: assignment(=) => assign_tabularEnergy
    generic   :: getInterF     => getInterF_withBin, getInterF_withoutBin

    procedure, private :: init_withPDF
    procedure, private :: init_withCDF
    procedure, private :: init_fromACE
    procedure, private :: assign_tabularEnergy
    procedure, private :: getInterF_withBin
    procedure, private :: getInterF_withoutBin
  end type tabularEnergy

contains

  !!
  !! Sample outgoing energy given Random Number Generator
  !!
  function sample(self, rand, bin) result (E)
    class(tabularEnergy), intent(in)         :: self
    class(RNG), intent(inout)                :: rand
    integer(shortInt), intent(out), optional :: bin
    real(defReal)                            :: E
    real(defReal)                            :: r

    r = rand % get()
    E = self % pdf % sample(r, bin)

  end function sample

  !!
  !! Return energy bounds of the probability distribution
  !!
  subroutine bounds(self,E_min,E_max)
    class(tabularEnergy), intent(in) :: self
    real(defReal), intent(out)       :: E_min
    real(defReal), intent(out)       :: E_max

    call self % pdf % bounds(E_min, E_max)

  end subroutine bounds

  !!
  !! Return probability that neutron was emitted with energy E
  !!
  function probabilityOf(self, E, bin) result(prob)
    class(tabularEnergy), intent(in)         :: self
    real(defReal), intent(in)                :: E
    integer(shortInt), intent(out), optional :: bin
    real(defReal)                            :: prob

    prob = self % pdf % probabilityOf(E)

  end function probabilityOf

  !!
  !! Initialise tabularEnergy with PDF only. Calculate CDF with PDF
  !!
  subroutine init_withPDF(self,E,PDF,interFlag)
    class(tabularEnergy), intent(inout)   :: self
    real(defReal),dimension(:),intent(in) :: E,PDF
    integer(shortInt),intent(in)          :: interFlag
    character(100),parameter              :: Here='init (tabularEnergy_class.f90)'

    if(count( E < 0.0 ) > 0) call fatalError(Here,'E contains -ve values')

    call self % pdf % init(E,PDF,interFlag)

  end subroutine init_withPDF

  !!
  !! Initialise tabularEnergy with PDF and CDF
  !! Does NOT check if PDF and CDF are consistant
  !!
  subroutine init_withCDF(self,E,PDF,CDF,interFlag)
    class(tabularEnergy), intent(inout)   :: self
    real(defReal),dimension(:),intent(in) :: E, PDF, CDF
    integer(shortInt),intent(in)          :: interFlag
    character(100),parameter              :: Here='init (tabularEnergy_class.f90)'


    if(count( E < 0.0 ) > 0) call fatalError(Here,'E contains -ve values')

    call self % pdf % init(E,PDF,CDF,interFlag)

  end subroutine init_withCDF

  !!
  !! Initialise tabularEnergy from ACE
  !! Head of aceCard needs to be set to the beginning of energy pdf data
  !! Uses the CDF in ACE data to initialise.
  !!
  subroutine init_fromACE(self, ACE)
    class(tabularEnergy), intent(inout)     :: self
    class(aceCard), intent(inout)           :: ACE
    integer(shortInt)                       :: INTT
    integer(shortInt)                       :: NP
    real(defReal),dimension(:), allocatable :: eGrid
    real(defReal),dimension(:), allocatable :: pdf
    real(defReal),dimension(:), allocatable :: cdf
    character(100),parameter :: Here = 'init_fromACE (tabularEnergy_class.f90)'

    ! Read data from ACE card
    INTT = ACE % readInt()

    ! Call error if interpolation flag indicates photon lines
    if( INTT > 10 ) then
      call fatalError(Here,'INTT > 10. Discrete photons lines are not yet implemented')

    end if

    ! Read rest of the data
    NP    = ACE % readInt()         ! Number of points
    eGrid = ACE % readRealArray(NP) ! Energy Values
    pdf   = ACE % readRealArray(NP) ! Probability density function
    cdf   = ACE % readRealArray(NP) ! Cumulative distribution function

    ! Initialise
    call self % init(eGrid, pdf, cdf, INTT)

  end subroutine init_fromACE

  !!
  !! Assignment
  !!
  subroutine assign_tabularEnergy(LHS,RHS)
    class(tabularEnergy),intent(out) :: LHS
    type(tabularEnergy),intent(in)   :: RHS

    LHS % pdf = RHS % pdf

  end subroutine assign_tabularEnergy

  !!
  !! Return interpolation factor for energy value E
  !! If E is outside bounds behaviour is undefined
  !!
  elemental function getInterF_withoutBin(self, E) result(f)
    class(tabularEnergy), intent(in) :: self
    real(defReal), intent(in)        :: E
    real(defReal)                    :: f

    f = self % pdf % getInterF(E)

  end function getInterF_withoutBin

  !!
  !! Return interpolation factor for energy value E, with bin number provided
  !! DOES NOT check if the bin is correct
  !!
  elemental function getInterF_withBin(self, E, bin) result(f)
    class(tabularEnergy), intent(in) :: self
    real(defReal), intent(in)        :: E
    integer(shortInt), intent(in)    :: bin
    real(defReal)                    :: f

    f = self % pdf % getInterF(E, bin)

  end function getInterF_withBin

  !!
  !! Constructor from PDF only
  !!
  function new_tabularEnergy(E,PDF,interFlag) result (new)
    real(defReal),dimension(:),intent(in) :: E,PDF
    integer(shortInt),intent(in)          :: interFlag
    type(tabularEnergy)                   :: new

    call new % init(E,PDF,interFlag)

  end function new_tabularEnergy


  !!
  !! Constructor form PDF and CDF
  !!
  function new_tabularEnergy_withCDF(E,PDF,CDF,interFlag) result (new)
    real(defReal),dimension(:),intent(in) :: E, PDF, CDF
    integer(shortInt),intent(in)          :: interFlag
    type(tabularEnergy)                   :: new

    call new % init(E,PDF,CDF,interFlag)

  end function new_tabularEnergy_withCDF

  !!
  !! Constructor from ACE
  !! Head of aceCard needs to be set to the beginning of energy pdf data
  !! Uses the CDF in ACE data to initialise.
  !!
  function new_tabularEnergy_fromACE(ACE) result (new)
    class(aceCard), intent(inout)            :: ACE
    type(tabularEnergy)                     :: new

    ! Initialise
    call new % init_fromACE(ACE)

  end function new_tabularEnergy_fromACE

end module tabularEnergy_class
