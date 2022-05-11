module thermalScatteringData_class

  use numPrecision
  use endfConstants
  use genericProcedures,            only : fatalError, numToChar, binarySearch, &
                                           endfInterpolate, isSorted
  use dataDeck_inter,               only : dataDeck
  use RNG_class,                    only : RNG
  use aceSabCard_class,             only : aceSabCard

  use uncorrelatedReactionCE_inter,  only : uncorrelatedReactionCE
  use thermalScatterInelastic_class, only : thInelasticScatter
  use thermalScatterElastic_class,   only : thElasticScatter

  implicit none
  private

  !! Private type to store thermal scattering cross sections
  type, private :: scatteringTable
    real(defReal), dimension(:), allocatable :: eGrid
    real(defReal), dimension(:), allocatable :: xs
  end type scatteringTable

  !!
  !! Public type for cross sections and energy/angle distributions
  !!
  !! Reads the tables from the nuclide's ACE card.
  !!
  !! Public Members:
  !!   eGrid  -> Energy grid
  !!
  !!
  type, public :: thermalData
    class(uncorrelatedReactionCE), allocatable :: inelasticOut
    class(uncorrelatedReactionCE), allocatable :: elasticOut
    type(scatteringTable)  :: inelastic
    type(scatteringTable)  :: elastic
    logical(defBool)       :: hasElastic = .false.
    logical(defBool)       :: isCoherent

  contains

    procedure :: init
    procedure :: kill
    procedure :: getEbounds
    procedure :: getInelXS
    procedure :: getElXS
    procedure :: buildFromACE

  end type thermalData

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
    class(thermalData), intent(inout) :: self
    class(dataDeck), intent(inout)    :: data
    character(100), parameter :: Here = 'init (thermalScatteringTables_class.f90)'

    ! Select buld procedure approperiate for given dataDeck
    select type(data)
    type is (aceSabCard)
        call self % buildFromACE(data)

      class default
        call fatalError(Here,'Thermal scattering data cannot be build from '//data % myType())
    end select

    allocate(thInelasticScatter :: self % inelasticOut)
    call self % inelasticOut % init(data, N_N_ThermINEL)
    if (self % hasElastic) then
      allocate(thElasticScatter :: self % elasticOut)
      call self % elasticOut % init(data, N_N_ThermEL)
    end if

  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(thermalData), intent(inout) :: self

  end subroutine kill

  !!
  !! Return the boundary values of the energy grid
  !!
  function getEbounds(self, request) result(eBounds)
    class(thermalData), intent(in) :: self
    character(*), intent(in)       :: request
    real(defReal), dimension(2)    :: eBounds
    real(defReal)                  :: E1, E2
    integer(shortInt)              :: N1, N2
    character(100), parameter :: Here = 'getEbounds (thermalScatteringTables_class.f90)'

    select case(request)
      case('inelastic')

        N1 = size(self % inelastic % eGrid)
        E1 = self % inelastic % eGrid(1)
        E2 = self % inelastic % eGrid(N1)

      case('elastic')

        if (self % hasElastic) then
          E1 = 1.0E-11_defReal
          N1 = size(self % elastic % eGrid)
          N2 = size(self % inelastic % eGrid)
          E2 = max(self % elastic % eGrid(N1),self % inelastic % eGrid(N2))
        else
          E1 = ZERO
          E2 = ZERO
        end if

      case default
        call fatalError(Here,'Unrecognised request string: |' // request //'|')

    end select

    eBounds = (/E1,E2/)

  end function getEbounds

  !!
  !! Calculate inelastic cross section
  !!
  !! Args:
  !!   E [in]    -> Energy of ingoing neutron
  !!   xi [in]   -> Random number
  !!   val [out] -> Interpolated values from tables
  !!
  pure subroutine getInelXS(self, E, val)
    class(thermalData), intent(in) :: self
    real(defReal), intent(in)      :: E
    real(defReal), intent(out)     :: val
    real(defReal)                  :: f
    integer(shortInt)              :: idx

    ! Get energy indexes
    idx = binarySearch(self % inelastic % eGrid, E)

    associate(E_top => self % inelastic % eGrid(idx + 1), &
              E_low  => self % inelastic % eGrid(idx))
      f = (E - E_low) / (E_top - E_low)
    end associate

    val = self % inelastic % xs(idx + 1) * f + (ONE-f) * self % inelastic % xs(idx)

  end subroutine getInelXS

  !!
  !! Calculate elastic cross section
  !!
  !! Args:
  !!   E [in]    -> Energy of ingoing neutron
  !!   xi [in]   -> Random number
  !!   val [out] -> Interpolated values from tables
  !!
  pure subroutine getElXS(self, E, val)
    class(thermalData), intent(in) :: self
    real(defReal), intent(in)      :: E
    real(defReal), intent(out)     :: val
    real(defReal)                  :: f
    integer(shortInt)              :: idx, N

    if (self % hasElastic) then
      N = size(self % elastic % eGrid)

      if (E < self % elastic % eGrid(1)) then
        val = ZERO
      elseif (E > self % elastic % eGrid(N)) then
        if (self % isCoherent) then
          val = self % elastic % xs(N)/E
        else
          val = ZERO
        end if
      else
        ! Get energy indexes
        idx = binarySearch(self % elastic % eGrid, E)
        if (self % isCoherent) then
          val = self % elastic % xs(idx)/E
        else
          associate(E_top => self % elastic % eGrid(idx + 1), &
                    E_low  => self % elastic % eGrid(idx))
            f = (E - E_low) / (E_top - E_low)
          end associate
          val = self % elastic % xs(idx + 1)*f + (ONE-f) * self % elastic % xs(idx)
        end if
      end if

    else
      val = ZERO
    end if

  end subroutine getElXS

  !!
  !! Build thermalScatteringData from ACE dataCard
  !!
  !! If the CDF is not sorted, the CDF doesn't end with one or there are negative
  !! cross sections, the tables are switched off for that nuclide
  !!
  !! Args:
  !!   ACE [inout] -> ACE card
  !!
  subroutine buildFromACE(self, ACE)
    class(thermalData), intent(inout) :: self
    type(aceSabCard), intent(inout)   :: ACE

    ! Get energy grids and cross sections
    self % inelastic % eGrid = ACE % ESZ_inelastic('energyGrid')
    self % inelastic % xs = ACE % ESZ_inelastic('inelasticXS')

    self % hasElastic = ACE % hasElastic()

    if (self % hasElastic) then
      self % elastic % eGrid = ACE % ESZ_elastic('energyGrid')
      self % elastic % xs = ACE % ESZ_elastic('Pvalues')
      self % isCoherent = ACE % isCoherent()
    end if

  end subroutine buildFromACE

end module thermalScatteringData_class
