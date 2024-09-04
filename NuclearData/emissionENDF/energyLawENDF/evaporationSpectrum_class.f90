module evaporationSpectrum_class

  use numPrecision
  use genericProcedures,   only : fatalError, isSorted
  use aceCard_class,       only : aceCard
  use RNG_class,           only : RNG
  use endfTable_class,     only : endfTable
  use energyLawENDF_inter, only : energyLawENDF

  implicit none
  private

  !!
  !! Constructors interface
  !!
  interface evaporationSpectrum
    module procedure new_evaporationSpectrum
    module procedure new_evaporationSpectrum_Inter
    module procedure new_evaporationSpectrum_fromACE
  end interface

  !!
  !! Outgoing energy spectrum given by EvaporationSpectrum (ACE LAW 9; ENDF LAW 9)
  !!   PDF(E_out) = C * E_out * exp(-E_out/T(E_in) )
  !!
  !!   T(E_in) is tabulated with values in ACE with 0 <= E_out <= E_in - U
  !!   U -> Restriction energy
  !!
  !! In order to sample this dictribution we make substitution E' = E_out/T(E_in) to obtain
  !!  PDF(E') = C' * E' exp(-E') = Gamma(2,1)
  !!
  !! We know we can sample Gamma(2,1) by noting it is equal to the distribution of sum
  !! of two Gamma(1,1) (Exponential) distributions.
  !!  PDF(E') = Gamma(1,1) + Gamma(1,1)
  !!
  !! Let X be a sample of PDF(E') and U & V uniformly distributed random numbers on [0,1]
  !!  X = -ln(U) - ln(V) = - ln(U*V)
  !!
  !! So sampled E_out = -T(E_in) * ln(U*V).
  !!
  !! It is also necessary to apply restriction in maximum energy, so E_out needs to be resampled if
  !!  E_out > E_in - U
  !!
  !!
  type, public, extends(energyLawENDF) :: evaporationSpectrum
    private
    type(endfTable)           :: T_of_E     ! "Temperature" as function of incident energy
    real(defReal)             :: U          ! Restriction energy
  contains
    ! Superclass procedures
    procedure :: sample
    procedure :: probabilityOf
    procedure :: kill

    ! Local procedures
    procedure :: init

  end type evaporationSpectrum

contains

  !!
  !! Sample outgoing energy given random number generator and incedent energy
  !!
  function sample(self,E_in,rand) result (E_out)
    class(evaporationSpectrum), intent(in) :: self
    real(defReal), intent(in)              :: E_in
    class(RNG), intent(inout)              :: rand
    real(defReal)                          :: E_out
    real(defReal)                          :: T, r1, r2

    ! Get value of T at energy E_in
    T = self % T_of_E % at(E_in)

    rejection:do
      r1 = rand % get()
      r2 = rand % get()

      E_out = -T * log(r1*r2)
      if(E_out <= E_in - self % U) exit rejection

    end do rejection

  end function sample

  !!
  !! Give probability of outgoing energy given incedent energy
  !!
  function probabilityOf(self,E_out,E_in) result (prob)
    class(evaporationSpectrum), intent(in) :: self
    real(defReal), intent(in)              :: E_out
    real(defReal), intent(in)              :: E_in
    real(defReal)                          :: prob
    real(defReal)                          :: T, C, U, Ep

    ! Get value of T at energy E_in
    T = self % T_of_E % at(E_in)

    ! Calculate normalisation constant
    U = self % U
    Ep = (E_in - U) /T
    C = T * T * (ONE - exp(-Ep) * (1+Ep))
    C = ONE /C

    ! Calculate probability
    prob = C * E_out * exp(-E_out/T)

  end function probabilityOf

  !!
  !! Initialise
  !!
  subroutine init(self,eGrid,T,U,bounds,interENDF)
    class(evaporationSpectrum), intent(inout)          :: self
    real(defReal),dimension(:),intent(in)              :: eGrid     ! T energy grid
    real(defReal),dimension(:),intent(in)              :: T         ! T values
    real(defReal), intent(in)                          :: U         ! Restriction energy
    integer(shortInt),dimension(:),intent(in),optional :: bounds    ! Bounds of interpolation regions
    integer(shortInt),dimension(:),intent(in),optional :: interENDF ! Interpolation flag
    character(100),parameter              :: Here='init (evaporationSpectrum_class.f90)'

    ! Perform sanity checks
    if(size(eGrid) /= size(T))      call fatalError(Here,'eGrid and T have diffrent size')
    if(.not.(isSorted(eGrid)))      call fatalError(Here,'eGrid is not sorted ascending')
    if ( any( eGrid < 0.0 )  )      call fatalError(Here,'eGrid contains -ve values')

    if ( any( T < 0.0 ) )           call fatalError(Here,'Neutron temperature T has -ve values')

    ! Initialise
    if (present(bounds) .and. present(interENDF)) then
      call self % T_of_E % init(eGrid,T,bounds,interENDF)

    else if ( present(bounds) .or. present(interENDF)) then
      call fatalError(Here,'Either "bounds" or "interENDF" is not given')

    else
      call self % T_of_E % init(eGrid,T) ! Default single region lin-lin interpolation

    end if

    self % U = U

  end subroutine init

  !!
  !! Release memory
  !!
  elemental subroutine kill(self)
    class(evaporationSpectrum), intent(inout) :: self

    call self % T_of_E % kill()
    self % U = ZERO

  end subroutine kill


  !!
  !! Simple constructor
  !! Only one lin-lin interpolation region
  !!
  function new_evaporationSpectrum(eGrid, T, U) result(new)
    real(defReal),dimension(:),intent(in)   :: eGrid
    real(defReal),dimension(:),intent(in)   :: T
    real(defReal), intent(in)               :: U
    type(evaporationSpectrum)                   :: new

    call new % init(eGrid,T,U)

  end function new_evaporationSpectrum

  !!
  !! Constructor
  !! Multiple interpolation regions & interpolation schemes
  !!
  function new_evaporationSpectrum_Inter(eGrid, T, U, bounds, interENDF) result(new)
    real(defReal),dimension(:),intent(in)     :: eGrid
    real(defReal),dimension(:),intent(in)     :: T
    real(defReal), intent(in)                 :: U
    integer(shortInt),dimension(:),intent(in) :: bounds, interENDF
    type(evaporationSpectrum)                     :: new

    call new % init(eGrid,T,U,bounds,interENDF)

  end function new_evaporationSpectrum_Inter

  !!
  !! Constructor from ACE
  !! aceCard head needs to be set to the beginning of the data
  !! Note: Should be accompanied by another initialisation routine to avoid reallocation of memory
  !!
  function new_evaporationSpectrum_fromACE(ACE) result(new)
    type(aceCard), intent(inout)               :: ACE
    type(evaporationSpectrum)                  :: new
    real(defReal),dimension(:),allocatable     :: eGrid
    real(defReal),dimension(:),allocatable     :: T
    real(defReal)                              :: U
    integer(shortInt),dimension(:),allocatable :: bounds
    integer(shortInt),dimension(:),allocatable :: interENDF
    integer(shortInt)                          :: NR
    integer(shortInt)                          :: N
    logical(defBool)                           :: hasInterRegions

    ! Read number of interpolation regions
    NR = ACE % readInt()
    hasInterRegions = (NR /= 0)

    ! Read interpolation data if present
    if ( hasInterRegions ) then
      bounds    = ACE % readIntArray(NR)
      interENDF = ACE % readIntArray(NR)

    end if

    ! Read rest of the data
    N     = ACE % readInt()        ! Number of incident energy points
    eGrid = ACE % readRealArray(N) ! Read incident energy grid
    T     = ACE % readRealArray(N) ! Temperature values
    U     = ACE % readReal()       ! Read restriction energy

    ! Initialise
    if (hasInterRegions) then
      call new % init(eGrid,T,U,bounds,interENDF)

    else
      call new % init(eGrid,T,U)

    end if

  end function new_evaporationSpectrum_fromACE



end module evaporationSpectrum_class
