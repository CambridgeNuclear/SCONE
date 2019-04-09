module multipleEnergyLaws_class
  use numPrecision
  use genericProcedures,       only : fatalError, numToChar
  use aceCard_class,           only : aceCard
  use RNG_class,               only : RNG
  use energyLawENDF_inter,     only : energyLawENDF
  use energyLawENDFSlot_class, only : energyLawENDFSlot
  use endfTable_class,         only : endfTable


  implicit none
  private

  !!
  !! Local helper type to bound ENDF Table with its bounds
  !!
  !! Public members:
  !!   table -> ENDF table
  !!   E_min -> Lower energy bounds
  !!   E_max -> Upper energy bounds
  !!
  type, private :: tableWrap
    type(endfTable) :: table
    real(defReal)   :: E_min
    real(defReal)   :: E_max
  end type tableWrap

  !!
  !! ENDF Energy Law to wrap multiple laws into one
  !!
  !! Sometimes a given reaction is represented as a mixture of multiple ENDF energy laws.
  !! Then each law is assigned with an energy dependent probability. Thus when sampling
  !! outgoing energy we need to roll a random number do decide which ENDF law to use.
  !! Each energy dependent probability[p(E)] is stored in a standard ENDF-style Table.
  !! If energy is below or above the bounds of the table, value at minimum or maximum is to
  !! be used.
  !! Specification taken from MCNP Manual Appendix F TABLE F-14 b)
  !!
  !!
  !! Private members:
  !!   prob -> Array of types that store ENDF p(E) table and its bounds
  !!   laws -> Array of energyLawENDF slots to store polymorphic instances of individual
  !!       energy laws.
  !!   num -> Number of already loaded energy laws
  !!
  !! Interface:
  !!   sample -> get outgoing energy given incident energy and random number
  !!   probabilityOf -> given incident and outgoing energy returns probability
  !!
  type, public,extends(energyLawENDF) :: multipleEnergyLaws
    private
    type(tableWrap), dimension(:), allocatable         :: prob
    type(energyLawENDFSlot), dimension(:), allocatable :: laws
    integer(shortInt)                                  :: num = 0
  contains
    ! Interface procedures
    procedure :: sample
    procedure :: probabilityOf

    ! Instance procedures
    procedure :: init
    procedure :: addLaw

  end type multipleEnergyLaws

contains


  !!
  !! Sample outgoing energy
  !!
  !! Args:
  !!   E_in [in] -> incident energy [MeV]
  !!   rand [inout] -> random number generator
  !!
  !! Returns:
  !!   Outgoing energy [MeV]
  !!
  !! Errors:
  !!   Gives fatalError if sampling fails to sample energy Law
  !!
  function sample(self, E_in, rand) result (E_out)
    class(multipleEnergyLaws), intent(in) :: self
    real(defReal), intent(in)             :: E_in
    class(RNG), intent(inout)             :: rand
    real(defReal)                         :: E_out
    real(defReal)                         :: r, E, prob
    integer(shortInt)                     :: i
    character(100), parameter :: Here = 'sample (multipleEnergyLaws_class.f90)'

    ! Generate random number
    r = rand % get()

    ! Find law index and sample
    do i =1,self % num
      ! Get probibility
      E = E_in
      E = max(E, self % prob(i) % E_min)
      E = min(E, self % prob(i) % E_max)
      prob = self % prob(i) % table % at(E)

      ! Check acceptance probability
      if ( r < prob) then
        E_out = self % laws(i) % sample(E_in, rand)
        return
      end if
    end do

    ! Error message
    call fatalError(Here, 'Failed to sample an energy law')
    E_out = ZERO

  end function sample

  !!
  !! Return probability of transition
  !!
  !! Args:
  !!   E_out [in] -> outgoing energy energy [MeV]
  !!   E_in [in] -> incident energy [MeV]
  !!
  !! Returns:
  !!   Probability of transition in [0,1]
  !!
  function probabilityOf(self,E_out,E_in) result (prob)
    class(multipleEnergyLaws), intent(in) :: self
    real(defReal), intent(in)             :: E_out,E_in
    real(defReal)                         :: prob
    real(defReal)                         :: E
    integer(shortInt)                     :: i

    prob = ZERO

    ! Calculate probability
    do i =1,self % num
      E = E_in
      E = max(E, self % prob(i) % E_min)
      E = min(E, self % prob(i) % E_max)

      prob = prob + self % laws(i) % probabilityOf(E_out, E_in) * self % prob(i) % table % at(E)

    end do

  end function probabilityOf

  !!
  !! Initialise multipleEnergyLaws
  !!
  !! Reserve space for N laws
  !!
  !! Args:
  !!   N [in] -> integer number of laws to reserve space for
  !!
  !! Errors:
  !!   FatalError if N <= 0
  !!
  subroutine init(self, N)
    class(multipleEnergyLaws), intent(inout) :: self
    integer(shortInt), intent(in)            :: N
    character(100),parameter :: Here ='init (multipleEnergyLaws_class.f90)'

    ! Verify N
    if ( N <= 0) then
      call fatalError(Here, 'Invalid number of energy Laws. Must be +ve. Is: ' // numToChar(N))
    end if

    ! Deallocate if allocated
    if(allocated(self % prob)) deallocate(self % prob)
    if(allocated(self % laws )) deallocate(self % laws)

    ! Allocate space
    allocate(self % prob(N))
    allocate(self % laws(N))

    ! Set number of loaded laws to 0
    self % num = 0

  end subroutine init

  !!
  !! Add additional energy law
  !!
  !! Add an energy law given already build polymorphic law and the data
  !! that specify its associated probability table.
  !!
  !! Args:
  !!   law [inout] -> allocated polymorphic energy law. Will be deallocated on exit
  !!   eGrid [in]  -> energy grid for probability table
  !!   pdf [in]    -> array of probability values associated with eGeid
  !!   bounds [in] -> optional integer array of bounds of diffrent ENDF
  !!       interpolation regions
  !!   interENDF [in] -> optional integer array of ENDF interpolation flags
  !!
  !! Errors:
  !!   Gives fatal error if:
  !!       * pdf or eGrid contains -ve entries
  !!       * law is not allocated
  !!       * one of bounds and interENDF is given without the other
  !!       * Maximum number of laws was already loaded
  !!
  subroutine addLaw(self, law, eGrid, pdf, bounds, interENDF)
    class(multipleEnergyLaws), intent(inout)             :: self
    class(energyLawENDF),allocatable, intent(inout)      :: law
    real(defReal), dimension(:), intent(in)              :: eGrid
    real(defReal), dimension(:), intent(in)              :: pdf
    integer(shortInt), dimension(:), intent(in),optional :: bounds
    integer(shortInt), dimension(:), intent(in),optional :: interENDF
    character(100), parameter :: Here = 'addLaw (multipleEnergyLaws_class.f90)'

    ! Check if space is availible
    self % num = self % num + 1
    if( self % num > size(self % prob)) then
      call fatalError(Here, 'Cannot add another law. Maximum was already reached')
    end if

    ! Check for -ve entries
    if(any(eGrid < ZERO) .or. any(pdf < 0)) then
      call fatalError(Here, '-ve entries in eGrid or pdf values ')
    end if

    ! Check if energy law is allocated
    if(.not.allocated(law)) then
      call fatalError(Here, 'Unallocated energy law was given.')
    end if

    ! Verify if bounds and interENDF were provided
    ! Build endfTable
    if(present(bounds) .and. present(interENDF)) then
      call self % prob(self % num) % table % init(eGrid, pdf, bounds, interENDF)

    else if( present(bounds) .eqv. present(interENDF)) then
      call self % prob(self % num) % table % init(eGrid, pdf)

    else
      call fatalError(Here,'Either bounds or interENDF was given without the other')

    end if
    self % prob(self % num) % E_min = eGrid(1)
    self % prob(self % num) % E_max = eGrid(size(eGrid))

    ! Move energy Law to slot
    call self % laws(self % num) % moveAllocFrom(law)

  end subroutine addLaw



end module multipleEnergyLaws_class
