module testFilter_test

  use numPrecision
  use testFilter_class, only : testFilter
  use particle_class,     only : particleState
  use dictionary_class,   only : dictionary
  use funit

  implicit none

@testCase
  type, extends(TestCase) :: test_testFilter
    private
    type(testFilter) :: filter
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_testFilter


contains

  !!
  !! Sets up test_testFilter object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_testFilter), intent(inout) :: this
    type(dictionary)                      :: tempDict

    call tempDict % init(2)
    call tempDict % store('minIdx', 4)
    call tempDict % store('maxIdx', 6)

    call this % filter % init(tempDict)

  end subroutine setUp

  !!
  !! Kills test_testFilter object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_testFilter), intent(inout) :: this

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test correct behaviour of the filter
  !!
@Test
  subroutine testFiltering(this)
    class(test_testFilter), intent(inout) :: this
    type(particleState)                   :: state
    logical(defBool)                      :: filterRes

    state % matIdx = 1
    filterRes = this % filter % isPass(state)
    @assertFalse(filterRes)

    state % matIdx = 4
    filterRes = this % filter % isPass(state)
    @assertTrue(filterRes)

    state % matIdx = 6
    filterRes = this % filter % isPass(state)
    @assertTrue(filterRes)

    state % matIdx = 7
    filterRes = this % filter % isPass(state)
    @assertFalse(filterRes)

  end subroutine testFiltering

end module testFilter_test
