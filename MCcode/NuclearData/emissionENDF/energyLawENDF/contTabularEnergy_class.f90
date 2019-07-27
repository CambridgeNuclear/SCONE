module contTabularEnergy_class

  use numPrecision
  use genericProcedures,   only : binarySearch, fatalError, interpolate, searchError, isSorted
  use aceCard_class,       only : aceCard
  use tabularEnergy_class, only : tabularEnergy
  use RNG_class,           only : RNG
  use energyLawENDF_inter, only : energyLawENDF


  implicit none
  private

  interface contTabularEnergy
    module procedure new_contTabularEnergy
    module procedure new_contTabularEnergy_fromACE
  end interface

  !!
  !! Continuous Tabular Distribution
  !! Energy distribution for an outgoing particle. LAW 4 in ENDF and ACE
  !! Does not support multiple interpolation regions & discrete photon lines
  !!
  type, public,extends(energyLawENDF):: contTabularEnergy
    private
    real(defReal),dimension(:),allocatable       :: eGrid
    type(tabularEnergy),dimension(:),allocatable :: ePdfs
  contains
    procedure :: sample
    procedure :: probabilityOf
    procedure :: kill

    procedure :: init

  end type contTabularEnergy

contains

  !!
  !! Sample outgoing particle energy given Random Number Generator and collision energy
  !!
  function sample(self,E_in,rand) result (E_out)
    class(contTabularEnergy), intent(in) :: self
    real(defReal), intent(in)            :: E_in
    class(RNG), intent(inout)            :: rand
    real(defReal)                        :: E_out
    integer(shortInt)                    :: idx
    real(defReal)                        :: r, eps
    real(defReal)                        :: E_min_low, E_max_low
    real(defReal)                        :: E_min_up, E_max_up
    real(defReal)                        :: E_min, E_max
    real(defReal)                        :: factor
    character(100),parameter             :: Here='sample (contTabularEnergy_class.f90)'

    idx = binarySearch(self % eGrid,E_in)
    call searchError(idx,Here)

    eps     = (E_in - self % eGrid(idx)) / (self % eGrid(idx+1) - self % eGrid(idx))

    r = rand % get()

    ! NOTE: Unlike in the equivalent tabular distribution for mu bounds of the
    !       energy distribution change between bins. We need to interpolate them
    !       as the result.

    ! Obtain bounds of the lower and upper bin
    call self % ePdfs(idx) % bounds(E_min_low, E_max_low)
    call self % ePdfs(idx) % bounds(E_min_up,  E_max_up )

    ! Calculate interpolated bounds
    E_min = E_min_low * (ONE - eps) + eps * E_min_up
    E_max = E_max_low * (ONE - eps) + eps * E_max_up

    ! Calculate interpolation between bounds of the distribution from which
    ! outgoing energy was sampled
    if(r < eps) then
      E_out = self % ePdfs(idx+1) % sample(rand)
      factor = (E_out- E_min_up)/(E_max_up - E_min_up)

    else
      E_out = self % ePdfs(idx) % sample(rand)
      factor = (E_out- E_min_low)/(E_max_low - E_min_low)

    end if

    ! Interpolate outgoing energy
    E_out = E_min *(ONE - factor) + factor * E_max

  end function sample

  !!
  !! Returns probability of emission of particle with energy E_out
  !! Given collision happend at energy E_in
  !! MAY NOT BE UP TO DATE
  function probabilityOf(self,E_out,E_in) result (prob)
    class(contTabularEnergy), intent(in) :: self
    real(defReal), intent(in)            :: E_out, E_in
    real(defReal)                        :: prob
    integer(shortInt)                    :: idx
    real(defReal)                        :: prob_1, prob_0, E_1, E_0
    character(100),parameter             :: Here='probabilityOf (contTabularEnergy_class.f90)'

    idx = binarySearch(self % eGrid,E_in)
    call searchError(idx,Here)

    prob_0 = self % ePdfs(idx)   % probabilityOf(E_out)
    prob_1 = self % ePdfs(idx+1) % probabilityOf(E_out)

    E_0 = self % eGrid(idx)
    E_1 = self % eGrid(idx+1)

    prob = interpolate(E_0, E_1, prob_0, prob_1, E_in)

  end function probabilityOf

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(contTabularEnergy), intent(inout) :: self

    if(allocated(self % eGrid)) deallocate(self % eGrid)

    if(allocated(self % ePdfs)) then
      call self % ePdfs % kill()
      deallocate(self % ePdfs)
    end if

  end subroutine kill


  !!
  !! Initialise from energy grid and table of single energy probabilty distributions
  !!
  subroutine init(self,eGrid,ePdfs)
    class(contTabularEnergy), intent(inout)     :: self
    real(defReal),dimension(:),intent(in)       :: eGrid
    type(tabularEnergy),dimension(:),intent(in) :: ePdfs
    character(100),parameter                    :: Here='init (contTabularEnergy_class.f90)'

    ! Check if the provided eGrid and ePdfs match in size and if eGrid is sorted and all its
    ! elements are +ve.
    if(size(eGrid) /= size(ePdfs))  call fatalError(Here,'eGrid and ePdfs have diffrent size')
    if(.not.(isSorted(eGrid)))      call fatalError(Here,'eGrid is not sorted ascending')
    if ( count( eGrid < 0.0 ) > 0 ) call fatalError(Here,'eGrid contains -ve values')


    if(allocated(self % eGrid)) deallocate(self % eGrid)
    if(allocated(self % ePdfs)) deallocate(self % ePdfs)

    self % eGrid  = eGrid
    self % ePdfs = ePdfs


  end subroutine init

  !!
  !! Constructor
  !!
  function new_contTabularEnergy(eGrid,ePdfs) result(new)
    real(defReal),dimension(:),intent(in)        :: eGrid
    type(tabularEnergy),dimension(:),intent(in)  :: ePdfs
    type(contTabularEnergy)                      :: new

    call new % init(eGrid,ePdfs)

  end function new_contTabularEnergy

  !!
  !! Constructor from ACE
  !! Head of aceCard needs to be set to the beginning of the data
  !! NOTE : Defining another init for ACE would help to avoid unnecesary reallocation of memory
  !!
  function new_contTabularEnergy_fromACE(ACE) result(new)
    type(aceCard), intent(inout)                 :: ACE
    type(contTabularEnergy)                      :: new
    integer(shortInt)                            :: NR
    integer(shortInt)                            :: N, i
    real(defReal),dimension(:),allocatable       :: eGrid
    integer(shortInt),dimension(:),allocatable   :: locEne
    type(tabularEnergy),dimension(:),allocatable :: ePdfs
    character(100),parameter :: Here = 'new_contTabularEnergy_fromACE (contTabularEnergy_class.f90)'

    ! Read number of interpolation regions
    NR = ACE % readInt()

    ! Return error if there are multiple regions
    if(NR /= 0) then
      call fatalError(Here,'Many inter. regions on energy distr. table are not supported')
    end if

    ! Read rest of the data
    N      = ACE % readInt()        ! Number of energy points
    eGrid  = ACE % readRealArray(N) ! Incident neutron energy grid
    locEne = ACE % readIntArray(N)  ! Read locations of PDF tables at a given incident energy

    allocate(ePdfs(N))

    ! Loop over all locations and read PDF at the given energy
    do i=1,N
      call ACE % setToEnergyLaw(locEne(i))
      ePdfs(i) = tabularEnergy(ACE)
    end do

    ! Initialise
    call new % init(eGrid,ePdfs)

  end function new_contTabularEnergy_fromACE


end module contTabularEnergy_class
