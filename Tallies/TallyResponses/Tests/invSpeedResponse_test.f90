module invSpeedResponse_test

  use numPrecision
  use universalVariables,     only : neutronMass, lightSpeed
  use invSpeedResponse_class, only : invSpeedResponse
  use particle_class,         only : particle, P_NEUTRON, P_PHOTON
  use dictionary_class,       only : dictionary
  use nuclearDatabase_inter,  only : nuclearDatabase
  use funit

  implicit none

@testCase
  type, extends(TestCase) :: test_invSpeedResponse
    private
    type(invSpeedResponse) :: response
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_invSpeedResponse


contains

  !!
  !! Sets up test_invSpeedResponse object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_invSpeedResponse), intent(inout) :: this
    type(dictionary)                            :: tempDict

  end subroutine setUp

  !!
  !! Kills test_invSpeedResponse object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_invSpeedResponse), intent(inout) :: this

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test correct behaviour of the response
  !!
@Test
  subroutine invSpeedResponseing(this)
    class(test_invSpeedResponse), intent(inout) :: this
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

  end subroutine invSpeedResponseing

end module invSpeedResponse_test
