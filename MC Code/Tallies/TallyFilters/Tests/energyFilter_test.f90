module energyFilter_test

  use numPrecision
  use energyFilter_class, only : energyFilter
  use particle_class,     only : particleState
  use dictionary_class,   only : dictionary
  use pFUnit_mod

  implicit none

@testCase
  type, extends(TestCase) :: test_energyFilter
    private
    type(energyFilter) :: filter
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_energyFilter

  !! Parameters
  real(defReal), parameter :: E_MIN = 1.27E-6_defReal
  real(defReal), parameter :: E_MAX = 1.34_defReal
  real(defReal), parameter :: E_Delta = 0.999_defReal ! Rel. diff. Must be < 1 for correct tests !


contains

  !!
  !! Sets up test_energyFilter object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_energyFilter), intent(inout) :: this
    type(dictionary)                        :: tempDict

    call tempDict % init(2)
    call tempDict % store('Emin', E_MIN)
    call tempDict % store('Emax', E_MAX)

    call this % filter % init(tempDict)

  end subroutine setUp

  !!
  !! Kills test_energyFilter object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_energyFilter), intent(inout) :: this

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test correct behaviour of the filter
  !!
@Test
  subroutine testFilter(this)
    class(test_energyFilter), intent(inout) :: this
    type(particleState)                     :: state
    real(defReal)                           :: testE
    logical(defBool)                        :: filterRes


    ! Below Emin -> FALSE
    testE = E_MIN * E_Delta
    state % E = testE
    filterRes = this % filter % filter(state)
    @assertFalse(filterRes)

    ! Emin -> TRUE
    testE = E_MIN
    state % E = testE
    filterRes = this % filter % filter(state)
    @assertTrue(filterRes)

    ! Below E max -> TRUE
    testE = E_MAX * E_Delta
    state % E = testE
    filterRes = this % filter % filter(state)
    @assertTrue(filterRes)

    ! E max -> TRUE
    testE = E_MAX
    state % E = testE
    filterRes = this % filter % filter(state)
    @assertTrue(filterRes)

    ! Above Emax -> FALSE
    testE = E_MAX / E_Delta
    state % E = testE
    filterRes = this % filter % filter(state)
    @assertFalse(filterRes)

  end subroutine testFilter

end module energyFilter_test
