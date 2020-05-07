module energyGrid_class

  use numPrecision
  use universalVariables
  use genericProcedures, only : fatalError, isDescending, binarySearch

  implicit none
  private

  !! Local parameters
  integer(shortInt), parameter :: UNDEF    = 0, &
                                  LIN      = 1, &
                                  LOGAR    = 2, &
                                  UNSTRUCT = 3
  !!
  !! Specialised type of energy grid for storing multi-group energy structures
  !!  Unlike grid_class it stores bin boundaries in descending order, Bin 1 has highest energy.
  !!  Interface:
  !!    bin(<int>)     -> returns value of the top value of energy group. Value at bin(NG+1) is
  !!                      the bottom of the lowest energy group (NG stands for number of groups)
  !!    search(<real>) -> returns the group number for the given value. Returns "valueOutsideArray"
  !!                      for a failed search.
  !!    getSize()      -> returns number of BINS(groups). Note that this is = [size(self % bins)-1].
  !!                      Returns 0 if grid is uninitialised
  !!
  !!  NOTE: internaly bins values are stored in reverse order (smallest at the begining)
  !!
  type, public :: energyGrid
    private
    integer(shortInt)                      :: type = UNDEF
    real(defReal),dimension(:),allocatable :: bins
    real(defReal)                          :: step = ZERO
  contains
    ! Initialisation procedures
    generic   :: init => init_equalSpaced, init_unstruct
    procedure :: init_equalSpaced
    procedure :: init_unstruct
    procedure :: kill

    ! Access procedures
    procedure :: bin
    procedure :: search
    procedure :: getSize
  end type energyGrid

contains

  !!
  !! Initialise equaly linearly spaced energyGrid
  !!
  subroutine init_equalSpaced(self, mini, maxi, N,type)
    class(energyGrid), intent(inout) :: self
    real(defReal), intent(in)        :: mini
    real(defReal), intent(in)        :: maxi
    integer(shortInt), intent(in)    :: N
    character(*), intent(in)         :: type
    integer(shortInt)                :: i
    character(100), parameter :: Here = 'init_equalSpaced ( energyGrid_class.f90)'

    ! Verify input
    if ( N < 1) call fatalError(Here,'Number of bins must be +ve')
    if ((maxi-mini)/maxi < FP_REL_TOL)  then
      call fatalError(Here,'Minimum value must be smaller then maximum above realtive FP tolerance')
    end if
    if(any([mini,maxi] < ZERO)) call fatalError(Here,'Energy grid requested contains -ve energies')

    ! Allocate space for bin boundaries
    allocate(self % bins(N+1))

    select case(type)
      case('lin')
        ! Calculate bin width
        self % step = (maxi - mini)/N

        ! Assign bin boundaries
        self % bins(1) = mini
        do i= 2,N+1
          self % bins(i) = mini + (i-1) * self % step
        end do

        ! Select energyGrid type
        self % type = LIN

      case('log')
        ! Check boundaries
        if (mini <= 0) call fatalError(Here,'For logarithmic energyGrid minimum must be +ve')

        ! Calculate bin width
        self % step = log(maxi/mini)/N

        ! Assign bin boundaries
        self % bins(1) = mini
        do i=2,N+1
          self % bins(i) = self % bins(i-1) * exp(self % step)
        end do

        ! Select energyGrid type
        self % type = LOGAR

      case default
        call fatalError(Here,'energyGrid type must be lin or log. Not: '//type)
    end select
  end subroutine init_equalSpaced

  !!
  !! Initialise unstructured energyGrid
  !!
  subroutine init_unstruct(self,bins)
    class(energyGrid), intent(inout)       :: self
    real(defReal),dimension(:), intent(in) :: bins
    character(100), parameter              :: Here = 'init_unstruct ( energyGrid_class.f90)'

    ! Verify input
    if( .not.isDescending(bins)) call fatalError(Here,'Provided energyGrid is not sorted descending')
    if( size(bins) < 2) call fatalError(Here,'Empty array or array of size 1 was provided')
    if(any(bins < ZERO)) call fatalError(Here,'Energy grid requested contains -ve energies')

    ! Initialise
    self % bins = bins(size(bins):1:-1)
    self % type = UNSTRUCT

  end subroutine init_unstruct

  !!
  !! Kill grid. Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(energyGrid), intent(inout) :: self

    if(allocated(self % bins)) deallocate(self % bins)
    self % step = ZERO
    self % type = UNDEF

  end subroutine kill


  !!
  !! Return value of the bins uder idx (upper bin limit)
  !! RETURNS -huge(defRal) if idx is outside bounds
  !!
  elemental function bin(self,idx)
    class(energyGrid), intent(in)       :: self
    integer(shortInt), intent(in)       :: idx
    real(defReal)                       :: bin

    if (idx > 0 .and. idx <= size(self % bins)) then
      bin = self % bins(size(self % bins) - idx +1)
    else
      bin = -huge(ONE)
    end if

  end function bin

  !!
  !! Returns index of the value in th energyGrid
  !! Returns  valueOutsideArray if target is outside bounds
  !! Returns  valueOutsideArray if energyGrid is uninitialised
  !!
  !! NOTE: Behaviour for values very close to bin boundary with exception of the first element
  !!       of self % bins is UNDEFINED, due to numerical precision.
  !!       When creating structured energyGrid there is going to be some accumulation
  !!       of numerical error and boundaries are not exactly  where they should be.
  !!       Thus exact value of the bin boundary is not easy to predict.
  !!
  elemental function search(self, value) result(idx)
    class(energyGrid), intent(in)  :: self
    real(defReal),intent(in) :: value
    integer(shortInt)        :: idx

    idx = 0
    select case(self % type)
      case(LIN)
        idx = floor((value-self % bins(1))/self % step)+1

      case(LOGAR)
        idx = floor(log(value/self % bins(1))/self % step) + 1

      case(UNSTRUCT)
        idx = binarySearch(self % bins, value)

    end select

    ! Flip index becouse the self % bins is inverted
    idx = size(self % bins) - idx

    ! Check whether errors happend
    if(idx < 1 .or. idx >= size(self % bins)) idx = valueOutsideArray

  end function search

  !!
  !! Returns number of bins (groups) in the energyGrid
  !! Not that this is size(self % bins) - 1
  !! Returns 0 for uninitialised array
  !!
  elemental function getSize(self) result(s)
    class(energyGrid), intent(in) :: self
    integer(shortInt)             :: s

    if(allocated(self % bins)) then
      s = size(self % bins) -1
    else
      s = 0
    end if

  end function getSize

end module energyGrid_class
