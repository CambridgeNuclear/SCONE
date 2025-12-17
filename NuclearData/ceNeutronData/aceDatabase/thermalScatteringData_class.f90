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
  !! Public type for thermal scattering cross sections and energy/angle distributions
  !! Reads the tables from the nuclide's S(a,b) ACE card.
  !! NOTE: when S(a,b) tables are on for inelastic scattering only, the elastic scattering
  !!       cross section must be set to zero
  !!
  !! Public Members:
  !!   inelasticOut -> reaction handle for inelastic scattering
  !!   elasticOut   -> reaction handle for elastic scattering
  !!   inelastic    -> tables for inelastic scattering
  !!   elastic      -> tables for elastic scattering
  !!   hasElastic   -> flag that indicates if elastic scattering is on
  !!   isCoherent   -> flag that indicates if elastic scatter is coherent or incoherent
  !!   kT           -> evaluation temperature in MeV
  !!
  !! Class Procedures:
  !!   init           -> initialises scattering tables
  !!   kill           -> returns to uninitialised state
  !!   getEbounds     -> returns energy grid boundaries for a required reaction
  !!   getInelXS      -> returns inelastic scattering xs
  !!   getElXS        -> returns elastic scattering xs
  !!   getTemperature -> returns the evaluation temperature
  !!   buildFromACE   -> build data from ACE card
  !!
  type, public :: thermalData
    class(uncorrelatedReactionCE), allocatable :: inelasticOut
    class(uncorrelatedReactionCE), allocatable :: elasticOut
    type(scatteringTable)  :: inelastic
    type(scatteringTable)  :: elastic
    logical(defBool)       :: hasElastic = .false.
    logical(defBool)       :: isCoherent = .false.
    real(defReal)          :: kT = -ONE

  contains

    procedure :: init
    procedure :: kill
    procedure :: getEbounds
    procedure :: getInelXS
    procedure :: getElXS
    procedure :: getTemperature
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

    ! initalise inelastic scattering
    allocate(thInelasticScatter :: self % inelasticOut)
    call self % inelasticOut % init(data, N_N_ThermINEL)

    ! initialise elastic scattering if data is present
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

    self % hasElastic = .false.
    self % isCoherent = .false.
    self % kT = -ONE

    call self % inelasticOut % kill()
    call self % elasticOut % kill()

    deallocate(self % inelastic % eGrid)
    deallocate(self % inelastic % xs)
    if(allocated(self % elastic % eGrid)) deallocate(self % elastic % eGrid)
    if(allocated(self % elastic % xs)) deallocate(self % elastic % xs)

  end subroutine kill

  !!
  !! Return the boundary values of the energy grid. Called only during initialisation
  !!
  !! Args:
  !!   request [in] -> string characterising which scattering reaction is required
  !!
  !! Errors:
  !! fatalerror if string other than 'inelastic' or 'elastic' is used
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
          ! Needed because if S(a,b) inelastic scattering data is present and
          ! elastic data is not, the elastic xs must be set to zero
          E2 = max(self % elastic % eGrid(N1),self % inelastic % eGrid(N2))
        else
          E1 = ZERO
          E2 = ZERO
        end if

      case default
        call fatalError(Here,'Unrecognised request string: |' // request //'|')

    end select

    eBounds = [E1, E2]

  end function getEbounds

  !!
  !! Return the temperature of the evaluated S(a,b) data
  !! Returns as kT - thermal energy in MeV
  !!
  pure function getTemperature(self) result(kT)
    class(thermalData), intent(in) :: self
    real(defReal)                  :: kT

    kT = self % kT 

  end function getTemperature

  !!
  !! Calculate inelastic cross section, returns cross section values
  !!
  !! Args:
  !!   E [in]    -> Energy of ingoing neutron
  !!
  pure function getInelXS(self, E) result(val)
    class(thermalData), intent(in) :: self
    real(defReal), intent(in)      :: E
    real(defReal)                  :: val
    real(defReal)                  :: f
    integer(shortInt)              :: idx

    ! Get energy indexes
    idx = binarySearch(self % inelastic % eGrid, E)

    associate(E_top => self % inelastic % eGrid(idx + 1), &
              E_low  => self % inelastic % eGrid(idx))
      f = (E - E_low) / (E_top - E_low)
    end associate

    val = self % inelastic % xs(idx + 1) * f + (ONE - f) * self % inelastic % xs(idx)

  end function getInelXS

  !!
  !! Calculate elastic cross section, returns cross section values
  !! It's called when S(a,b) tables are switched on, even if there's no thermal
  !! elastic scattering
  !!
  !! Args:
  !!   E [in]    -> Energy of ingoing neutron
  !!
  pure function getElXS(self, E) result(val)
    class(thermalData), intent(in) :: self
    real(defReal), intent(in)      :: E
    real(defReal)                  :: val
    real(defReal)                  :: f
    integer(shortInt)              :: idx, N

    ! Cross section is set to zero if the elastic scattering flag is false
    ! and if energy values is lower than the nuclide energy grid first value
    if (.not. self % hasElastic) then
      val = ZERO
      return
    elseif (E < self % elastic % eGrid(1)) then
      val = ZERO
      return
    end if

    ! Retrieve energy grid size
    N = size(self % elastic % eGrid)

    ! Coherent elastic scattering
    if (self % isCoherent) then

      if (E > self % elastic % eGrid(N)) then
        val = self % elastic % xs(N)/E
      else
        idx = binarySearch(self % elastic % eGrid, E)
        val = self % elastic % xs(idx)/E
      end if

    else  ! Incoherent elastic scattering

      if (E > self % elastic % eGrid(N)) then
        val = ZERO
      else
        idx = binarySearch(self % elastic % eGrid, E)
        associate(E_top => self % elastic % eGrid(idx + 1), &
                  E_low  => self % elastic % eGrid(idx))
          f = (E - E_low) / (E_top - E_low)
        end associate
        val = self % elastic % xs(idx + 1) * f + (ONE - f) * self % elastic % xs(idx)
      end if

    end if

  end function getElXS

  !!
  !! Build thermalScatteringData from ACE dataCard
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
    
    ! Read temperature
    self % kT = ACE % TZ

    ! Check if elastic scattering data is present
    self % hasElastic = ACE % hasElastic()

    ! Initialise elastic scattering
    if (self % hasElastic) then
      self % elastic % eGrid = ACE % ESZ_elastic('energyGrid')
      self % elastic % xs = ACE % ESZ_elastic('Pvalues')
      self % isCoherent = ACE % isCoherent()
    end if

  end subroutine buildFromACE

end module thermalScatteringData_class
