module fluxResponse_test

  use numPrecision
  use fluxResponse_class,    only : fluxResponse
  use particle_class,        only : particle
  use dictionary_class,      only : dictionary
  use nuclearDatabase_inter, only : nuclearDatabase
  use pFUnit_mod

  implicit none

@testCase
  type, extends(TestCase) :: test_fluxResponse
    private
    type(fluxResponse) :: response
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_fluxResponse


contains

  !!
  !! Sets up test_fluxResponse object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_fluxResponse), intent(inout) :: this
    type(dictionary)                      :: tempDict

  end subroutine setUp

  !!
  !! Kills test_fluxResponse object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_fluxResponse), intent(inout) :: this

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test correct behaviour of the filter
  !!
@Test
  subroutine fluxResponseing(this)
    class(test_fluxResponse), intent(inout) :: this
    type(particle)                          :: p
    class(nuclearDatabase),pointer          :: xsData

    @assertEqual(ONE, this % response % get(p, xsData), 1.0E-9_defReal)

  end subroutine fluxResponseing

end module fluxResponse_test
