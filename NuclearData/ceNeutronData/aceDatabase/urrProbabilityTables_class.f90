module urrProbabilityTables_class

  use numPrecision
  use endfConstants
  use genericProcedures,            only : fatalError, numToChar, binarySearch, &
                                           endfInterpolate, isSorted
  use RNG_class,                    only : RNG
  use dataDeck_inter,               only : dataDeck
  use aceCard_class,                only : aceCard

  implicit none
  private

  !! Private type to store probability tables
  type, private :: urrTable
    real(defReal), dimension(:), allocatable :: CDF
    real(defReal), dimension(:), allocatable :: tot
    real(defReal), dimension(:), allocatable :: el
    real(defReal), dimension(:), allocatable :: fiss
    real(defReal), dimension(:), allocatable :: capt
  end type urrTable

  !!
  !! Public type for unresolved resonance probability tables
  !!
  !! Reads the tables from the nuclide's ACE card. For now it doesn't save the heating
  !! number.
  !!
  !! It saves the IOA, which is not used in the code
  !!
  !! Public Members:
  !!   nGrid  -> Dimension of the energy grid
  !!   nTable -> Dimension of probability table
  !!   INT    -> Interpolation flag
  !!   ILF    -> Inelatic competition flag (only used to check if inelatic scattering is zero)
  !!   IOA    -> Other absorption flag. NOTE: not used in the code
  !!   IFF    -> Multiplication factor flag
  !!   eGrid  -> Energy grid
  !!   table  -> Array of probability tables
  !!
  !!
  type, public :: urrProbabilityTables
    integer(shortInt) :: nGrid
    integer(shortInt) :: nTable
    integer(shortInt) :: INT
    integer(shortInt) :: ILF
    integer(shortInt) :: IOA
    integer(shortInt) :: IFF
    real(defReal), dimension(:), allocatable  :: eGrid
    type(urrTable), dimension(:), allocatable :: table

  contains

    procedure :: init
    procedure :: kill
    procedure :: getEbounds
    procedure :: getIFF
    procedure :: sampleXSs
    procedure :: buildFromACE

  end type urrProbabilityTables

contains

  !!
  !! Initialise
  !!
  !! Args:
  !!   data[inout] -> cross sections data deck
  !!
  !! Errors:
  !!   fatalError if data type is not ACE
  !!
  subroutine init(self, data)
    class(urrProbabilityTables), intent(inout) :: self
    class(dataDeck), intent(inout)             :: data
    character(100), parameter :: Here = 'init (urrProbabilityTables_class.f90)'

    ! Select buld procedure approperiate for given dataDeck
    select type(data)
      type is (aceCard)
        call self % buildFromACE(data)

      class default
        call fatalError(Here,'Elastic neutron scattering cannot be build from '//data % myType())
    end select

  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(urrProbabilityTables), intent(inout) :: self

    if (allocated(self % eGrid)) deallocate(self % eGrid)
    if (allocated(self % table)) deallocate(self % table)

  end subroutine kill

  !!
  !! Return the boundary values of the energy grid
  !!
  function getEbounds(self) result(eBounds)
    class(urrProbabilityTables), intent(in) :: self
    real(defReal), dimension(2)             :: eBounds

    eBounds(1) = self % eGrid(1)
    eBounds(2) = self % eGrid(self % nGrid)

  end function getEbounds

  !!
  !! Return the multiplication factor flag
  !!
  function getIFF(self) result(IFF)
    class(urrProbabilityTables), intent(in) :: self
    integer(shortInt)                       :: IFF

    IFF = self % IFF

  end function getIFF

  !!
  !! Samples from probability tables
  !!
  !! Args:
  !!   E [in]    -> Energy of ingoing neutron
  !!   xi [in]   -> Random number
  !!   val [out] -> Interpolated values from tables
  !!
  subroutine sampleXSs(self, E, xi, val)
    class(urrProbabilityTables), intent(in)    :: self
    real(defReal), intent(in)                  :: E
    real(defReal), intent(in)                  :: xi
    real(defReal), dimension(3), intent(out)   :: val
    real(defReal)                              :: E1, E2, res1, res2
    integer(shortInt)                          :: enIdx1, enIdx2, idx1, idx2
    real(defReal), dimension(self % nTable +1) :: tmpCDF

    ! Get energy indexes
    enIdx1 = binarySearch(self % eGrid, E)
    enIdx2 = enIdx1 + 1

    ! Get energy values
    E1 = self % eGrid(enIdx1)
    E2 = self % eGrid(enIdx2)

    ! Add zero to CDF and get indexes
    tmpCDF(1) = ZERO
    tmpCDF(2:self % nTable +1) = self % table(enIdx1) % CDF
    idx1 = binarySearch(tmpCDF, xi)
    tmpCDF(2:self % nTable +1) = self % table(enIdx2) % CDF
    idx2 = binarySearch(tmpCDF, xi)

    ! Interpolate to get elastic cross section
    res1 = self % table(enIdx1) % el(idx1)
    res2 = self % table(enIdx2) % el(idx2)
    val(1) = endfInterpolate(E1,E2,res1,res2,E,self % INT)

    ! Interpolate to get capture (n,gamma) cross section
    res1 = self % table(enIdx1) % capt(idx1)
    res2 = self % table(enIdx2) % capt(idx2)
    val(2) = endfInterpolate(E1,E2,res1,res2,E,self % INT)

    ! Interpolate to get fission cross section
    res1 = self % table(enIdx1) % fiss(idx1)
    res2 = self % table(enIdx2) % fiss(idx2)
    val(3) = endfInterpolate(E1,E2,res1,res2,E,self % INT)

  end subroutine sampleXSs

  !!
  !! Build urrProbabilityTables from ACE dataCard
  !!
  !! If the CDF is not sorted, the CDF doesn't end with one or there are negative
  !! cross sections, the tables are switched off for that nuclide
  !!
  !! Args:
  !!   ACE [inout] -> ACE card
  !!
  subroutine buildFromACE(self, ACE)
    class(urrProbabilityTables), intent(inout) :: self
    type(aceCard), intent(inout)               :: ACE
    integer(shortInt)                          :: i

    ! Set head in ACE card
    call ACE % setToProbTab()

    self % nGrid  = ACE % readInt()        ! Read size of the energy grid
    self % nTable = ACE % readInt()        ! Read size of probability table
    self % INT    = ACE % readInt()        ! Read interpolation factor
    self % ILF    = ACE % readInt()        ! Read inelastic flag
    self % IOA    = ACE % readInt()        ! Read other absorptions flag
    self % IFF    = ACE % readInt()        ! Read multiplication factor flag

    allocate(self % eGrid(self % nGrid))
    self % eGrid = ACE % readRealArray(self % nGrid)  ! Read energy grid

    allocate(self % table(self % nGrid))

    ! Read probability tables
    do i = 1,self % nGrid
      self % table(i) % CDF  = ACE % readRealArray(self % nTable)
      self % table(i) % tot  = ACE % readRealArray(self % nTable)
      self % table(i) % el   = ACE % readRealArray(self % nTable)
      self % table(i) % fiss = ACE % readRealArray(self % nTable)
      self % table(i) % capt = ACE % readRealArray(self % nTable)
      call ACE % advanceHead(self % nTable)
    end do

    ! Discard this table if values don't make sense
    do i = 1,self % nGrid
      if (.not.isSorted(self % table(i) % CDF)) then
        print '(A)', "Probability table discarded because CDF is not sorted"
        call self % kill()
        return
      elseif(self % table(i) % CDF(self % nTable) /= ONE) then
        print '(A)', "Probability table discarded because CDF does not end with 1.0 "
        call self % kill()
        return
      elseif (self % table(i) % CDF(self % nTable) /= ONE .or. &
              any(self % table(i) % tot < ZERO) .or. any(self % table(i) % el < ZERO) .or. &
              any(self % table(i) % fiss < ZERO) .or. any(self % table(i) % capt < ZERO) ) then
        print '(A)', "Probability table discarded because negative cross-sections present "
        call self % kill()
        return
      end if
    end do

  end subroutine buildFromACE

end module urrProbabilityTables_class
