module kalbach87_class

  use numPrecision
  use genericProcedures,       only : fatalError, binarySearch, searchError, interpolate, isSorted
  use aceCard_class,           only : aceCard
  use RNG_class,               only : RNG
  use kalbachPdf_class,        only : kalbachPdf
  use correlatedLawENDF_inter, only : correlatedLawENDF

  implicit none
  private

  interface kalbach87
    module procedure new_kalbach87
    module procedure new_kalbach87_fromACE
  end interface

  !!
  !! Kalbach-87 Formalism (ACE LAW 44. ENDF File 6 Law 1)
  !! Correlated mu-energy distribution
  !! Does not support multiple interpolation regions on incident energy grid.
  !! Provisionaly tested.
  !! probabilityOf was NOT properly verified
  !!
  type, public,extends(correlatedLawENDF) :: kalbach87
    private
    real(defReal),dimension(:),allocatable    :: eGrid
    type(kalbachPdf),dimension(:),allocatable :: pdfs
  contains
    procedure :: sample
    procedure :: probabilityOf
    procedure :: kill

    procedure :: init

  end type kalbach87

contains

  !!
  !! Samples mu and E_out given incident energy E_in and random number generator
  !!
  subroutine sample(self,mu,E_out,E_in,rand)
    class(kalbach87), intent(in)  :: self
    real(defReal), intent(out)    :: mu
    real(defReal), intent(out)    :: E_out
    real(defReal), intent(in)     :: E_in
    class(RNG), intent(inout)     :: rand
    integer(shortInt)             :: idx
    real(defReal)                 :: E_min_low, E_max_low
    real(defReal)                 :: E_min_up, E_max_up
    real(defReal)                 :: E_min, E_max
    real(defReal)                 :: factor
    real(defReal)                 :: r, eps
    character(100),parameter      :: Here='sample (kalbach87_class.f90)'

    ! Find Interval index
    idx = binarySearch(self % eGrid,E_in)
    call searchError(idx,Here)

    ! Calculate threshold and sample random number
    eps = (E_in - self % eGrid(idx)) / (self % eGrid(idx+1) - self % eGrid(idx))

    r = rand % get()

    ! NOTE: Unlike in the equivalent tabular distribution for mu bounds of the
    !       energy distribution change between bins. We need to interpolate them
    !       as the result.

    ! Obtain bounds of the lower and upper bin
    call self % pdfs(idx)   % bounds(E_min_low, E_max_low)
    call self % pdfs(idx+1) % bounds(E_min_up,  E_max_up )

    ! Calculate interpolated bounds
    E_min = E_min_low * (ONE - eps) + eps * E_min_up
    E_max = E_max_low * (ONE - eps) + eps * E_max_up

    ! Calculate interpolation between bounds of the distribution from which
    ! outgoing energy was sampled
    if(r < eps) then
      call self % pdfs(idx+1) % sample(mu, E_out, rand)
      factor = (E_out- E_min_up)/(E_max_up - E_min_up)

    else
      call self % pdfs(idx) % sample(mu, E_out, rand)
      factor = (E_out- E_min_low)/(E_max_low - E_min_low)

    end if
    
    ! Interpolate outgoing energy
    E_out = E_min *(ONE - factor) + factor * E_max

  end subroutine sample

  !!
  !! Returns probability that neutron was emitted at angle mu and with E_out
  !! given incident energy E_in.
  !!
  function probabilityOf(self,mu,E_out,E_in) result(prob)
    class(kalbach87), intent(in) :: self
    real(defReal), intent(in)    :: mu
    real(defReal), intent(in)    :: E_out
    real(defReal), intent(in)    :: E_in
    real(defReal)                :: prob
    integer(shortInt)            :: idx
    real(defReal)                :: prob_1, prob_0, E_1, E_0
    character(100),parameter     :: Here='probabilityOf (kalbach87_class.f90)'

    ! Find interval index
    idx = binarySearch(self % eGrid,E_in)
    call searchError(idx,Here)

    ! Obtain probabilities & energies at boundaries of the interval
    prob_0 = self % pdfs(idx)   % probabilityOf(mu,E_out)
    prob_1 = self % pdfs(idx+1) % probabilityOf(mu,E_out)

    E_0 = self % eGrid(idx)
    E_1 = self % eGrid(idx+1)

    ! Interpolate
    prob = interpolate(E_0, E_1, prob_0, prob_1, E_in)

  end function probabilityOf

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(Self)
    class(kalbach87), intent(inout) :: self

    if(allocated(self % eGrid)) deallocate(self % eGrid)

    if(allocated(self % pdfs)) then
      call self % pdfs % kill()
      deallocate(self % pdfs)
    end if

  end subroutine kill


  !!
  !! Initialise
  !!
  subroutine init(self,eGrid,pdfs)
    class(kalbach87), intent(inout) :: self
    real(defReal),dimension(:)      :: eGrid
    type(kalbachPdf),dimension(:)   :: pdfs
    character(100),parameter        :: Here='init (kalbach87_class.f90)'

    ! Check if the provided eGrid and pdfs match in size and if eGrid is sorted and all its
    ! elements are +ve.
    if(size(eGrid) /= size(pdfs))   call fatalError(Here,'eGrid and ePdfs have diffrent size')
    if(.not.(isSorted(eGrid)))      call fatalError(Here,'eGrid is not sorted ascending')
    if(any( eGrid < 0.0 ))          call fatalError(Here,'eGrid contains -ve values')

    ! Deallocate current contents if allocated
    if(allocated(self % eGrid)) deallocate(self % eGrid)
    if(allocated(self % pdfs))  deallocate(self % pdfs)

    ! Assign content
    self % eGrid = eGrid
    self % pdfs  = pdfs

  end subroutine init

  !!
  !! Constructor
  !!
  function new_kalbach87(eGrid,pdfs) result(new)
    real(defReal),dimension(:),intent(in)    :: eGrid
    type(kalbachPdf),dimension(:),intent(in) :: pdfs
    type(kalbach87)                          :: new

    call new % init(eGrid,pdfs)

  end function new_kalbach87

  !!
  !! Constructor from ACE
  !! aceCard read head needs to be set to the beginning of the data
  !!
  function new_kalbach87_fromACE(ACE) result(new)
    type(aceCard), intent(inout)               :: ACE
    type(kalbach87)                            :: new
    real(defReal),dimension(:),allocatable     :: eGrid
    type(kalbachPdf),dimension(:),allocatable  :: pdfs
    integer(shortInt),dimension(:),allocatable :: locKal
    integer(shortInt)                          :: NR, N, i
    character(100),parameter :: Here ='new_kalbach87_fromACE (kalbach87_class.f90)'

    ! Read number of interpolation regions
    NR = ACE % readInt()

    ! Return error if there are multiple regions
    if(NR /= 0) then
      call fatalError(Here,'Many inter. regions on energy distr. table are not supported')
    end if

    ! Read rest of the data
    N      = ACE % readInt()        ! Number of energy points
    eGrid  = ACE % readRealArray(N) ! Incident neutron energy grid
    locKal = ACE % readIntArray(N)  ! Read locations of PDF tables at a given incident energy

    allocate(pdfs(N))

    ! Loop over all locations and read PDF at the given energy
    do i=1,N
      call ACE % setToEnergyLaw(locKal(i))
      pdfs(i) = kalbachPdf(ACE)
    end do

    ! Initialise
    call new % init(eGrid,pdfs)

  end function new_kalbach87_fromACE

end module kalbach87_class
