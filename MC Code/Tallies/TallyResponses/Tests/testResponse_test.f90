module testResponse_test

  use numPrecision
  use testResponse_class, only : testResponse
  use particle_class,     only : particle
  use dictionary_class,   only : dictionary
  use pFUnit_mod

  implicit none

@testCase
  type, extends(TestCase) :: test_testResponse
    private
    type(testResponse) :: response
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_testResponse


contains

  !!
  !! Sets up test_testResponse object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_testResponse), intent(inout) :: this
    type(dictionary)                        :: tempDict

    call tempDict % init(1)
    call tempDict % store('value', 1.3_defReal)
    call this % response % init(tempDict)

  end subroutine setUp

  !!
  !! Kills test_testResponse object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_testResponse), intent(inout) :: this

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test correct behaviour of the filter
  !!
@Test
  subroutine testResponseing(this)
    class(test_testResponse), intent(inout) :: this
    type(particle)                          :: p

    @assertEqual(1.3_defReal, this % response % get(p), 1.0E-9_defReal)

  end subroutine testResponseing

end module testResponse_test
