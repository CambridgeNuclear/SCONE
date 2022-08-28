module contTabularEnergy_class

  use numPrecision
  use endfConstants
  use genericProcedures,   only : binarySearch, fatalError, interpolate, searchError, isSorted,&
                                  ceilingSearch => linearCeilingIdxOpen_shortInt, numToChar
  use aceCard_class,       only : aceCard
  use tabularEnergy_class, only : tabularEnergy
  use RNG_class,           only : RNG
  use energyLawENDF_inter, only : energyLawENDF


  implicit none
  private

  interface contTabularEnergy
    module procedure new_contTabularEnergy_fromACE
  end interface

  ! Local parameters
  ! Locations of bounds & flags in the rank 2 inter array with interpolation data
  integer(shortInt), parameter :: INT_BOUNDS = 1
  integer(shortInt), parameter :: INT_FLAGS  = 2

  !!
  !! Continuous Tabular Distribution
  !!
  !! Energy distribution for an outgoing particle. LAW 4 in ENDF and ACE
  !! Discrete photon lines have not been implemented yet.
  !! Only lin-lin and histogram interpolation is supported
  !!
  !! Private Members:
  !!   eGrid  -> Energy Grid [MeV]
  !!   ePdfs  -> Array of probability distributions associated with each energy point
  !!   inter  -> Interpolation parameters
  !!
  !! Interface:
  !!   energyLawENDF interface
  !!   init -> Initialise from components
  !!
  type, public,extends(energyLawENDF):: contTabularEnergy
    private
    real(defReal),dimension(:),allocatable         :: eGrid
    type(tabularEnergy),dimension(:),allocatable   :: ePdfs
    integer(shortInt), dimension(:,:), allocatable :: inter
  contains
    ! Superclass interface
    procedure :: sample
    procedure :: probabilityOf
    procedure :: kill

    ! Local procedures
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
    integer(shortInt)                    :: idx, flag, inter_idx
    real(defReal)                        :: r, eps
    real(defReal)                        :: E_min_low, E_max_low
    real(defReal)                        :: E_min_up, E_max_up
    real(defReal)                        :: E_min, E_max
    real(defReal)                        :: factor
    character(100),parameter             :: Here = 'sample (contTabularEnergy_class.f90)'

    idx = binarySearch(self % eGrid,E_in)
    call searchError(idx,Here)

    ! Get interpolation flag
    if (allocated(self % inter)) then
      inter_idx = ceilingSearch(self % inter(:,INT_BOUNDS), idx)
      if (inter_idx <= 0) then
        call fatalError(Here, 'Failed interpolation region search: '//numToChar(inter_idx))
      end if
      flag = self % inter(inter_idx, INT_FLAGS)

    else
      flag = linLinInterpolation

    end if

    ! Select correct interpolation
    select case (flag)
      case (linLinInterpolation)
        ! NOTE: Unlike in the equivalent tabular distribution for mu bounds of the
        !       energy distribution change between bins. We need to interpolate them
        !       as the result.

        ! Obtain bounds of the lower and upper bin
        call self % ePdfs(idx) % bounds(E_min_low, E_max_low)
        call self % ePdfs(idx + 1) % bounds(E_min_up,  E_max_up )

        eps = (E_in - self % eGrid(idx)) / (self % eGrid(idx+1) - self % eGrid(idx))

        ! Calculate interpolated bounds
        E_min = E_min_low * (ONE - eps) + eps * E_min_up
        E_max = E_max_low * (ONE - eps) + eps * E_max_up

        ! Calculate interpolation between bounds of the distribution from which
        ! outgoing energy was sampled
        r = rand % get()
        if(r < eps) then
          E_out = self % ePdfs(idx+1) % sample(rand)
          factor = (E_out- E_min_up)/(E_max_up - E_min_up)

        else
          E_out = self % ePdfs(idx) % sample(rand)
          factor = (E_out- E_min_low)/(E_max_low - E_min_low)

        end if

        ! Interpolate outgoing energy
        E_out = E_min *(ONE - factor) + factor * E_max

      case (histogramInterpolation)
        ! Sample from bottom bin
        E_out = self % ePdfs(idx) % sample(rand)

      case default
        call fatalError(Here, 'Unsupported interpolation flag: '//numToChar(flag))
        E_out = ZERO ! Make compiler happy

      end select
  end function sample

  !!
  !! Returns probability of emission of particle with energy E_out
  !!
  !! See energyLawENDF_inter for details
  !!
  !! This function has not been implemented. `probabilityOf` will be
  !! removed soon from all emissionENDF.
  !!
  function probabilityOf(self,E_out,E_in) result (prob)
    class(contTabularEnergy), intent(in) :: self
    real(defReal), intent(in)            :: E_out, E_in
    real(defReal)                        :: prob
    ! integer(shortInt)                    :: idx
    ! real(defReal)                        :: prob_1, prob_0, E_1, E_0
    character(100),parameter             :: Here = 'probabilityOf (contTabularEnergy_class.f90)'

    call fatalError(Here, 'probabilityOf has not been implemented.')
    prob = ZERO

    ! idx = binarySearch(self % eGrid, E_in)
    ! call searchError(idx,Here)
    !
    ! prob_0 = self % ePdfs(idx)   % probabilityOf(E_out)
    ! prob_1 = self % ePdfs(idx+1) % probabilityOf(E_out)
    !
    ! E_0 = self % eGrid(idx)
    ! E_1 = self % eGrid(idx+1)
    !
    ! prob = interpolate(E_0, E_1, prob_0, prob_1, E_in)

  end function probabilityOf

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(contTabularEnergy), intent(inout) :: self

    if(allocated(self % eGrid)) deallocate(self % eGrid)
    if(allocated(self % inter)) deallocate(self% inter)

    if(allocated(self % ePdfs)) then
      call self % ePdfs % kill()
      deallocate(self % ePdfs)
    end if

  end subroutine kill


  !!
  !! Initialise from energy grid and table of single energy probabilty distributions
  !!
  !! Args:
  !!   eGrid [in]     -> Energy grid [MeV]
  !!   ePdfs [in]     -> Array of initialised outgoing energy probability distributions
  !!   bounds [in]    -> Optional. Array of tops of diffrent interpolation regions
  !!   interENDF [in] -> Optional. Corresponding interpolation flags
  !!
  !! Errors:
  !!   FatalError if:
  !!     * Size mismatch between eGrid & ePdfs
  !!     * eGrid is not sorted or contains -ve entries
  !!     * bounds or interENDF is given without the other
  !!     * bounds does not specify interpolation for the entire table
  !!
  subroutine init(self, eGrid, ePdfs, bounds, interENDF)
    class(contTabularEnergy), intent(inout)               :: self
    real(defReal), dimension(:), intent(in)               :: eGrid
    type(tabularEnergy), dimension(:), intent(in)         :: ePdfs
    integer(shortInt), dimension(:), optional, intent(in) :: bounds
    integer(shortInt), dimension(:), optional, intent(in) :: interENDF
    character(100), parameter :: Here = 'init (contTabularEnergy_class.f90)'

    ! Check if the provided eGrid and ePdfs match in size and if eGrid is sorted and all its
    ! elements are +ve.
    if(size(eGrid) /= size(ePdfs))  call fatalError(Here,'eGrid and ePdfs have diffrent size')
    if(.not.(isSorted(eGrid)))      call fatalError(Here,'eGrid is not sorted ascending')
    if ( count( eGrid < 0.0 ) > 0 ) call fatalError(Here,'eGrid contains -ve values')


    if(allocated(self % eGrid)) deallocate(self % eGrid)
    if(allocated(self % ePdfs)) deallocate(self % ePdfs)

    self % eGrid = eGrid
    self % ePdfs = ePdfs

    ! Check interpolation data
    if (present(bounds) .neqv. present(interENDF)) then
      call fatalError(Here, 'Either bounds or interENDF is given without the other')

    else if (present(bounds)) then
      if (size(bounds) /= size(interENDF)) then
        call fatalError(Here, 'Size of bounds and  interENDF does not match')

      else if (.not.isSorted(bounds)) then
        call fatalError(Here, 'Array with interpolation region bounds is not sorted ascending')

      else if (any(bounds < 0)) then
        call fatalError(Here, 'Bounds contain -ve entries')

      else if (maxval(bounds) > size(eGrid)) then
        call fatalError(Here, 'Bounds contains values larger then size of energy grid')

      else if (bounds(size(bounds)) /= size(eGrid)) then
        call fatalError(Here, 'Incomplete interpolation scheme')

      end if

      ! Load the interpolation data
      allocate(self % inter(size(bounds), 2))
      self % inter(:, INT_BOUNDS) = bounds
      self % inter(:, INT_FLAGS)  = interENDF

    end if

  end subroutine init

  !!
  !! Constructor from ACE
  !!
  !! NOTE:
  !!   Defining another init for ACE would help to avoid unnecesary reallocation of memory
  !!
  !! Args:
  !!   ACE [inout] -> Ace card set to the beggining of the data
  !!   root [in]   -> Root location on XSS of ACE for the relative locators of single-energy
  !!     probability tables
  !!
  !! Result:
  !!   Initialised instance of contTabularEnergy
  !!
  function new_contTabularEnergy_fromACE(ACE, root) result(new)
    type(aceCard), intent(inout)                 :: ACE
    integer(shortInt), intent(in)                :: root
    type(contTabularEnergy)                      :: new
    integer(shortInt)                            :: NR
    integer(shortInt)                            :: N, i
    real(defReal),dimension(:),allocatable       :: eGrid
    integer(shortInt),dimension(:),allocatable   :: locEne
    type(tabularEnergy),dimension(:),allocatable :: ePdfs
    character(100),parameter :: Here = 'new_contTabularEnergy_fromACE (contTabularEnergy_class.f90)'

    ! Read Interpolation data
    NR = ACE % readInt()
    associate (bounds => ACE % readIntArray(NR), flags => ACE % readIntArray(NR))

      ! Read rest of the data
      N      = ACE % readInt()        ! Number of energy points
      eGrid  = ACE % readRealArray(N) ! Incident neutron energy grid
      locEne = ACE % readIntArray(N)  ! Read locations of PDF tables at a given incident energy

      allocate(ePdfs(N))

      ! Loop over all locations and read PDF at the given energy
      do i = 1, N
        call ACE % setRelativeTo(root, locEne(i))
        ePdfs(i) = tabularEnergy(ACE)
      end do

      ! Initialise
      if (NR == 0) then
        call new % init(eGrid,ePdfs)

      else if (NR > 0) then
        call new % init(eGrid, ePdfs, bounds, flags)

      else
        call fatalError(Here, '-ve number of interpolation regions: '//numToChar(NR))

      end if
    end associate

  end function new_contTabularEnergy_fromACE


end module contTabularEnergy_class
