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
    type(energyFilter) :: mgFilter
    type(energyFilter) :: cemgFilter
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_energyFilter

  !! Parameters
  real(defReal), parameter    :: E_MIN = 1.27E-6_defReal
  real(defReal), parameter    :: E_MAX = 1.34_defReal
  real(defReal), parameter    :: E_Delta = 0.999_defReal ! Rel. diff. Must be < 1 for correct tests !
  integer(shortInt),parameter :: G_TOP = 3
  integer(shortInt),parameter :: G_LOW = 7

contains

  !!
  !! Sets up test_energyFilter object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_energyFilter), intent(inout) :: this
    type(dictionary)                        :: tempDict

    ! Build CE filter
    call tempDict % init(4)
    call tempDict % store('Emin', E_MIN)
    call tempDict % store('Emax', E_MAX)
    call this % filter % init(tempDict)

    ! Build MG-CE filter
    call tempDict % store('Gtop', G_TOP)
    call tempDict % store('Glow', G_LOW)
    call this % cemgFilter % init(tempDict)

    ! Build MG filter
    call tempDict % kill()
    call tempDict % init(2)
    call tempDict % store('Gtop', G_TOP)
    call tempDict % store('Glow', G_LOW)
    call this % mgFilter % init(tempDict)

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
  !! Test correct behaviour of the CE filter
  !!
@Test
  subroutine testCEFilter(this)
    class(test_energyFilter), intent(inout) :: this
    type(particleState)                     :: state
    real(defReal)                           :: testE
    logical(defBool)                        :: filterRes

    ! Below Emin -> FALSE
    testE = E_MIN * E_Delta
    state % E = testE
    filterRes = this % filter % isPass(state)
    @assertFalse(filterRes)

    ! Emin -> TRUE
    testE = E_MIN
    state % E = testE
    filterRes = this % filter % isPass(state)
    @assertTrue(filterRes)

    ! Below E max -> TRUE
    testE = E_MAX * E_Delta
    state % E = testE
    filterRes = this % filter % isPass(state)
    @assertTrue(filterRes)

    ! E max -> TRUE
    testE = E_MAX
    state % E = testE
    filterRes = this % filter % isPass(state)
    @assertTrue(filterRes)

    ! Above Emax -> FALSE
    testE = E_MAX / E_Delta
    state % E = testE
    filterRes = this % filter % isPass(state)
    @assertFalse(filterRes)

  end subroutine testCEFilter

  !!
  !! Test correct behaviour of the CEMG filter
  !!
@Test
  subroutine testMGFilter(this)
    class(test_energyFilter), intent(inout) :: this
    type(particleState)                     :: state
    integer(shortInt)                       :: testG
    logical(defBool)                        :: filterRes

    state % isMG = .true.
    ! Below G_TOP -> FALE
    testG = G_TOP - 1
    state % G = testG
    filterRes = this % cemgFilter % isPass(state)
    @assertFalse(filterRes)

    ! Exactly G_TOP -> TRUE
    testG = G_TOP
    state % G = testG
    filterRes = this % cemgFilter % isPass(state)
    @assertTrue(filterRes)

    ! Above G_TOP -> TRUE
    testG = G_TOP + 1
    state % G = testG
    filterRes = this % cemgFilter % isPass(state)
    @assertTrue(filterRes)

    ! Above G_LOW -> FALSE
    testG = G_LOW + 1
    state % G = testG
    filterRes = this % cemgFilter % isPass(state)
    @assertFalse(filterRes)

    ! Exactly G_LOW -> TRUE
    testG = G_LOW
    state % G = testG
    filterRes = this % cemgFilter % isPass(state)
    @assertTRUE(filterRes)

  end subroutine testMGFilter

  !!
  !! Test correct behaviour of the CEMG filter
  !!
@Test
  subroutine testCEMGFilter(this)
    class(test_energyFilter), intent(inout) :: this
    type(particleState)                     :: state
    real(defReal)                           :: testE
    integer(shortInt)                       :: testG
    logical(defBool)                        :: filterRes

    state % isMG = .false.
    ! Below Emin -> FALSE
    testE = E_MIN * E_Delta
    state % E = testE
    filterRes = this % cemgFilter % isPass(state)
    @assertFalse(filterRes)

    ! Emin -> TRUE
    testE = E_MIN
    state % E = testE
    filterRes = this % cemgFilter % isPass(state)
    @assertTrue(filterRes)

    ! Below E max -> TRUE
    testE = E_MAX * E_Delta
    state % E = testE
    filterRes = this % cemgFilter % isPass(state)
    @assertTrue(filterRes)

    ! E max -> TRUE
    testE = E_MAX
    state % E = testE
    filterRes = this % cemgFilter % isPass(state)
    @assertTrue(filterRes)

    ! Above Emax -> FALSE
    testE = E_MAX / E_Delta
    state % E = testE
    filterRes = this % cemgFilter % isPass(state)
    @assertFalse(filterRes)

    state % isMG = .true.
    ! Below G_TOP -> FALE
    testG = G_TOP - 1
    state % G = testG
    filterRes = this % cemgFilter % isPass(state)
    @assertFalse(filterRes)

    ! Exactly G_TOP -> TRUE
    testG = G_TOP
    state % G = testG
    filterRes = this % cemgFilter % isPass(state)
    @assertTrue(filterRes)

    ! Above G_TOP -> TRUE
    testG = G_TOP + 1
    state % G = testG
    filterRes = this % cemgFilter % isPass(state)
    @assertTrue(filterRes)

    ! Above G_LOW -> FALSE
    testG = G_LOW + 1
    state % G = testG
    filterRes = this % cemgFilter % isPass(state)
    @assertFalse(filterRes)

    ! Exactly G_LOW -> TRUE
    testG = G_LOW
    state % G = testG
    filterRes = this % cemgFilter % isPass(state)
    @assertTRUE(filterRes)

  end subroutine testCEMGFilter



end module energyFilter_test
