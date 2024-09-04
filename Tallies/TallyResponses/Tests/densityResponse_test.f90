module densityResponse_test

  use numPrecision
  use universalVariables,    only : neutronMass, lightSpeed
  use densityResponse_class, only : densityResponse
  use particle_class,        only : particle, P_NEUTRON, P_PHOTON
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

    ! Test neutron density with different particle energies
    p % type = P_NEUTRON
    p % isMG = .false.
    p % E = ONE
    res   = ONE / lightSpeed / sqrt(TWO * p % E / neutronMass)
    @assertEqual(res, this % response % get(p, xsData), res*1.0E-9_defReal)

    p % E = 1.6e-06_defReal
    res   = ONE / lightSpeed / sqrt(TWO * p % E / neutronMass)
    @assertEqual(res, this % response % get(p, xsData), res*1.0E-9_defReal)

    ! Test photon density
    p % type = P_PHOTON
    res   = ONE/lightSpeed
    @assertEqual(res, this % response % get(p, xsData), res*1.0E-9_defReal)

  end subroutine densityResponseing

end module densityResponse_test
