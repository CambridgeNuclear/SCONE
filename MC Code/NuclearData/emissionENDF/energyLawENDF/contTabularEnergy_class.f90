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
    character(100),parameter             :: Here='sample (contTabularEnergy_class.f90)'

    idx = binarySearch(self % eGrid,E_in)
    call searchError(idx,Here)

    eps = (E_in - self % eGrid(idx)) / (self % eGrid(idx+1) - self % eGrid(idx))
    r = rand % get()

    if(r < eps) then
      E_out = self % ePdfs(idx+1) % sample(rand)
    else
      E_out = self % ePdfs(idx) % sample(rand)
    end if

  end function sample

  !!
  !! Returns probability of emission of particle with energy E_out
  !! Given collision happend at energy E_in
  !!
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
