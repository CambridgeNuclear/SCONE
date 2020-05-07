module maxwellSpectrum_class

  use numPrecision
  use genericProcedures,      only : fatalError, linearFloorIdxClosed_Real, searchError, isSorted, &
                                     interpolate
  use aceCard_class,          only : aceCard
  use RNG_class,              only : RNG
  use maxwellEnergyPdf_class, only : maxwellEnergyPdf
  use energyLawENDF_inter,    only : energyLawENDF
  use endfTable_class,        only : endfTable

  implicit none
  private

  integer(shortInt), parameter :: maxIter = 1000

  interface searchGrid
    module procedure linearFloorIdxClosed_Real
  end interface

  interface maxwellSpectrum
    module procedure new_maxwellSpectrum
    module procedure new_maxwellSpectrum_Inter
    module procedure new_maxwellSpectrum_fromACE
  end interface

  !!
  !! Simple Maxwell Fission Spectrum (ENDF LAW 7 )
  !!
  type, public, extends(energyLawENDF) :: maxwellSpectrum
    private
    type(maxwellEnergyPdf)    :: maxwellPdf  ! Maxwell Energy pdf
    type(endfTable)           :: T_of_E      ! "Temperature" as function of incident energy
    real(defReal)             :: U   = ZERO  ! Restriction energy

  contains
    procedure :: sample
    procedure :: probabilityOf
    procedure :: kill

    procedure :: init

  end type maxwellSpectrum

contains

  !!
  !! Samples energy of outgoing neutron given incedent particle energy and random number generator
  !!
  function sample(self,E_in,rand) result (E_out)
    class(maxwellSpectrum), intent(in) :: self
    real(defReal),intent(in)           :: E_in
    class(RNG),intent(inout)           :: rand
    real(defReal)                      :: E_out
    integer(shortInt)                  :: i
    real(defReal)                      :: T
    character(100),parameter           :: Here='sample (maxwellSpectrum_class.f90)'

    ! Interpolate value of T
    T = self % T_of_E % at(E_in)

    ! Sample Maxwellian distribution
    do i=1,maxIter

      E_out = self % maxwellPdf % sample(T,rand)
      if (E_out < E_in - self% U) return

    end do

    ! Write error messege if sampling was in infinate loop
    call fatalError(Here,'Sampling failed to be accepted after "maxIter" iterations')

  end function sample

  !!
  !! Returns probability that particle was emmited with energy E_out given incident energy E_in
  !!
  function probabilityOf(self,E_out,E_in) result (prob)
    class(maxwellSpectrum), intent(in) :: self
    real(defReal), intent(in)          :: E_out
    real(defReal), intent(in)          :: E_in
    real(defReal)                      :: prob
    real(defReal)                      :: T, C
    real(defReal)                      :: Arg, Sqrt_Arg, inv_C
    ! Interpolate T
    T = self % T_of_E % at(E_in)

    ! Calculate normalisation constant
    Arg      = (E_in-self % U) / T
    Sqrt_Arg = sqrt(Arg)

    inv_C = sqrt(T)*T * (0.5*sqrt(PI)*erf(Sqrt_Arg) - Sqrt_Arg * exp(-Arg))
    C = 1/inv_C

    ! Obtain probability
    prob = self % maxwellPdf % probabilityOf(E_out,T,C)

  end function probabilityOf

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(maxwellSpectrum), intent(inout) :: self

    call self % T_of_E % kill()
    self % U = ZERO

  end subroutine kill

  !!
  !! Initialise
  !!
  subroutine init(self,eGrid,T,U,bounds,interENDF)
    class(maxwellSpectrum), intent(inout)              :: self
    real(defReal),dimension(:),intent(in)              :: eGrid     ! T energy grid
    real(defReal),dimension(:),intent(in)              :: T         ! T values
    real(defReal), intent(in)                          :: U         ! Restriction energy
    integer(shortInt),dimension(:),intent(in),optional :: bounds    ! Bounds of interpolation regions
    integer(shortInt),dimension(:),intent(in),optional :: interENDF ! Interpolation flag
    character(100),parameter              :: Here='init (maxwellSpectrum_class.f90)'

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
  !! Simple constructor
  !! Only one lin-lin interpolation region
  !!
  function new_maxwellSpectrum(eGrid,T,U) result(new)
    real(defReal),dimension(:),intent(in)   :: eGrid
    real(defReal),dimension(:),intent(in)   :: T
    real(defReal), intent(in)               :: U
    type(maxwellSpectrum)                   :: new

    call new % init(eGrid,T,U)

  end function new_maxwellSPectrum
    
  !!
  !! Constructor
  !! Multiple interpolation regions & interpolation schemes
  !!
  function new_maxwellSpectrum_Inter(eGrid,T,U,bounds,interENDF) result(new)
    real(defReal),dimension(:),intent(in)     :: eGrid
    real(defReal),dimension(:),intent(in)     :: T
    real(defReal), intent(in)                 :: U
    integer(shortInt),dimension(:),intent(in) :: bounds, interENDF
    type(maxwellSpectrum)                     :: new

    call new % init(eGrid,T,U,bounds,interENDF)

  end function new_maxwellSPectrum_Inter

  !!
  !! Constructor from ACE
  !! aceCard head needs to be set to the beginning of the data
  !! Note: Should be accompanied by another initialisation routine to avoid reallocation of memory
  !!
  function new_maxwellSpectrum_fromACE(ACE) result(new)
    type(aceCard), intent(inout)               :: ACE
    type(maxwellSpectrum)                      :: new
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

  end function new_maxwellSpectrum_fromACE


end module maxwellSpectrum_class
