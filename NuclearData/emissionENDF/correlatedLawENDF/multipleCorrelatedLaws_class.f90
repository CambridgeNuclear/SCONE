module multipleCorrelatedLaws_class

  use numPrecision
  use genericProcedures,           only : fatalError, numToChar
  use aceCard_class,               only : aceCard
  use RNG_class,                   only : RNG
  use correlatedLawENDF_inter,     only : correlatedLawENDF
  use correlatedLawENDFslot_class, only : correlatedLawENDFslot
  use endfTable_class,             only : endfTable

  implicit none
  private

  !!
  !! Local helper type to bound ENDF Table with its bounds
  !!
  !! Public Members:
  !!   table -> ENDF table
  !!   E_min -> Lower energy bounds
  !!   E_max -> Upper energy bounds
  !!
  type, private :: tableWrap
    type(endfTable) :: table
    real(defReal)   :: E_min = ZERO
    real(defReal)   :: E_max = ZERO
  end type tableWrap

  !!
  !! ENDF Correlated Law to wrap multiple laws into one
  !!
  !! Motivation and implementation is based on equivalent `multipleEnergyLaws`. Please refer
  !! to itsdocumentation for extra details/references.
  !!
  !! Private Members:
  !!   prob -> Array of `tableWraps` to store energy-dependent probablity of a given law [p(E)]
  !!   laws -> Array of correlatedLawENDF slots to store polymorphic instances of individual laws.
  !!   num  -> Number of loaded correlated laws
  !!
  !! Interface:
  !!   correlatedLawENDF interface
  !!   init   -> Initialise saving space for N correlated laws
  !!   addLaw -> Add correlated law
  !!
  !! Usage:
  !!   To initialise
  !!
  type, public, extends(correlatedLawENDF) :: multipleCorrelatedLaws
    private
    type(tableWrap), dimension(:), allocatable             :: prob
    type(correlatedLawENDFslot), dimension(:), allocatable :: laws
    integer(shortInt)                                      :: num = 0
  contains
    ! correlatedLawENDF interface
    procedure :: sample
    procedure :: probabilityOf
    procedure :: kill

    ! Class procedures
    procedure :: init
    procedure :: addLaw

  end type multipleCorrelatedLaws

contains

  !!
  !! Sample outgoing energy & angle
  !!
  !! See correlatedLawENDF_inter for details
  !!
  subroutine sample(self, mu, E_out, E_in, rand)
    class(multipleCorrelatedLaws), intent(in) :: self
    real(defReal), intent(out)                :: mu
    real(defReal), intent(out)                :: E_out
    real(defReal), intent(in)                 :: E_in
    class(RNG), intent(inout)                 :: rand
    real(defReal)                             :: r, E, prob
    integer(shortInt)                         :: i
    character(100), parameter :: Here = 'sample (multipleCorrelatedLaws_class.f90)'

    ! Generate random number
    r = rand % get()

    ! Find law index and sample
    do i = 1, self % num
      ! Get probibility
      E = E_in
      E = max(E, self % prob(i) % E_min)
      E = min(E, self % prob(i) % E_max)
      prob = self % prob(i) % table % at(E)

      ! Check acceptance probability
      if ( r < prob) then
        call self % laws(i) % sample(mu, E_out, E_in, rand)
        return

      end if
      ! Decrement roll value
      r = r - prob
    end do

    ! Error message
    call fatalError(Here, 'Failed to sample correlated law. Check that total probability sums to 1?')
    mu    = -TWO ! Make compiler happy
    E_out = ZERO

  end subroutine sample

  !!
  !! Return probability of outgoing mu & E_out
  !!
  !! See correlatedLawENDF_inter for details
  !!
  function probabilityOf(self, mu, E_out, E_in) result(prob)
    class(multipleCorrelatedLaws), intent(in) :: self
    real(defReal), intent(in)                 :: mu
    real(defReal), intent(in)                 :: E_out
    real(defReal), intent(in)                 :: E_in
    real(defReal)                             :: prob
    real(defReal)                             :: E
    integer(shortInt)                         :: i

    prob = ZERO

    ! Calculate probability
    do i = 1, self % num
      E = E_in
      E = max(E, self % prob(i) % E_min)
      E = min(E, self % prob(i) % E_max)

      prob = prob + self % laws(i) % probabilityOf(mu, E_out, E_in) * self % prob(i) % table % at(E)

    end do

  end function probabilityOf

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(multipleCorrelatedLaws), intent(inout) :: self

    if(allocated(self % prob)) then
      call self % prob % table % kill()
      deallocate(self % prob)
    end if

    if(allocated(self % laws)) then
      call self % laws % kill()
      deallocate(self % laws)
    end if

  end subroutine kill

  !!
  !! Initialise multipleCorrelatedLaws
  !!
  !! Reserve space for N laws
  !!
  !! Args:
  !!   N [in] -> Integer. Number of laws that will be loaded into the object
  !!
  !! Errors:
  !!   FatalError if N <= 0
  !!
  subroutine init(self, N)
    class(multipleCorrelatedLaws), intent(inout) :: self
    integer(shortInt), intent(in)                :: N
    character(100), parameter :: Here = 'init (multipleCorrelatedLaws_class.f90)'

    ! Make sure is clean
    call self % kill()

    ! Verify N
    if (N <= 0) then
      call fatalError(Here, 'Invalid number of correlated laws. Must be +ve. Is: '//numToChar(N))
    end if

    ! Allocate space
    allocate(self % prob(N))
    allocate(self % laws(N))

    ! Set number of loaded laws
    self % num = 0

  end subroutine init

  !!
  !! Add additional correlated law
  !!
  !! Add a correlatedlaw and data for its associated p(E) table
  !!
  !! Args:
  !!   law [inout]    -> Allocated polymorphic instance of correlated law. Will be
  !!     DEALLOCATED on exit.
  !!   eGrid [in]     -> Energy grid for p(E) table
  !!   pdf [in]       -> Values of p(E) on the eGrid
  !!   bounds [in]    -> Optional. Array of bounds of diffrent ENDF interpolation regions.
  !!   interENDF [in] -> Optional. Interpolation flags assoctaed with the bounds
  !!
  !! Errors:
  !!   FatalError if:
  !!     * pdf or eGrid contain -ve entries or have diffrent length
  !!     * law is not allocated
  !!     * bounds & interENDF have diffrent length or one is missing
  !!     * Maximum number of laws was already loaded
  !!
  subroutine addLaw(self, law, eGrid, pdf, bounds, interENDF)
    class(multipleCorrelatedLaws), intent(inout)          :: self
    class(correlatedLawENDF), allocatable, intent(inout)  :: law
    real(defReal), dimension(:), intent(in)               :: eGrid
    real(defReal), dimension(:), intent(in)               :: pdf
    integer(shortInt), dimension(:), intent(in), optional :: bounds
    integer(shortInt), dimension(:), intent(in), optional :: interENDF
    character(100), parameter :: Here = 'addLaw (multipleCorrelatedLaws_class.f90)'

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



end module multipleCorrelatedLaws_class
