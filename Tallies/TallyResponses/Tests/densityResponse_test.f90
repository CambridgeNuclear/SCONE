module densityResponse_test

  use numPrecision
  use densityResponse_class, only : densityResponse
  use particle_class,        only : particle
  use dictionary_class,      only : dictionary
  use nuclearDatabase_inter, only : nuclearDatabase
  use pFUnit_mod

  implicit none

@testCase
  type, extends(TestCase) :: test_densityResponse
    private
    type(densityResponse) :: response
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_densityResponse


contains

  !!
  !! Sets up test_densityResponse object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_densityResponse), intent(inout) :: this
    type(dictionary)                           :: tempDict

  end subroutine setUp

  !!
  !! Kills test_densityResponse object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_densityResponse), intent(inout) :: this

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test correct behaviour of the response
  !!
@Test
  subroutine densityResponseing(this)
    class(test_densityResponse), intent(inout) :: this
    type(particle)                             :: p
    class(nuclearDatabase),pointer             :: xsData
    real(defReal)                              :: res

    ! Test with different particle energies
    p % E = ONE
    res   = ONE/1.3831592645e+09_defReal
    @assertEqual(res, this % response % get(p, xsData), res*1.0E-9_defReal)

    p % E = 1.6e-06_defReal
    res   = ONE/1.7495734571e+06_defReal
    @assertEqual(res, this % response % get(p, xsData), res*1.0E-9_defReal)

  end subroutine densityResponseing

end module densityResponse_test
