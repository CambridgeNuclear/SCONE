module macroResponse_test

  use numPrecision
  use endfConstants
  use macroResponse_class,            only : macroResponse
  use particle_class,                 only : particle, P_NEUTRON
  use dictionary_class,               only : dictionary
  use testNeutronDatabase_class,      only : testNeutronDatabase
  use funit

  implicit none

@testCase
  type, extends(TestCase) :: test_macroResponse
    private
    type(macroResponse)        :: response_total
    type(macroResponse)        :: response_capture
    type(macroResponse)        :: response_fission
    type(macroResponse)        :: response_nuFission
    type(macroResponse)        :: response_absorbtion
    type(macroResponse)        :: response_MT
    type(testNeutronDatabase)  :: xsData
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_macroResponse


contains

  !!
  !! Sets up test_macroResponse object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_macroResponse), intent(inout) :: this
    type(dictionary)                         :: tempDict

    ! Allocate and initialise test nuclearData
    ! Cross-sections:         Total         eScattering  IeScatter Capture      Fission      nuFission
    call this % xsData % build(6.0_defReal, 3.0_defReal, ZERO,     2.0_defReal, 1.0_defReal, 1.5_defReal)

    ! Set up responses
    ! Total
    call tempDict % init(2)
    call tempDict % store('type','macroResponse')
    call tempDict % store('MT', macroTotal)
    call this % response_total % init(tempDict)
    call tempDict % kill()

    ! Capture
    call tempDict % init(2)
    call tempDict % store('type','macroResponse')
    call tempDict % store('MT', macroCapture)
    call this % response_capture % init(tempDict)
    call tempDict % kill()

    ! Fission
    call tempDict % init(2)
    call tempDict % store('type','macroResponse')
    call tempDict % store('MT', macroFission)
    call this % response_fission % init(tempDict)
    call tempDict % kill()

    ! nuFission
    call tempDict % init(2)
    call tempDict % store('type','macroResponse')
    call tempDict % store('MT', macroNuFission)
    call this % response_nuFission % init(tempDict)
    call tempDict % kill()

    ! Absorbtion
    call tempDict % init(2)
    call tempDict % store('type','macroResponse')
    call tempDict % store('MT', macroAbsorbtion)
    call this % response_absorbtion % init(tempDict)
    call tempDict % kill()

    ! MT number
    call tempDict % init(2)
    call tempDict % store('type','macroResponse')
    call tempDict % store('MT', 105)
    call this % response_MT % init(tempDict)
    call tempDict % kill()

  end subroutine setUp

  !!
  !! Kills test_macroResponse object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_macroResponse), intent(inout) :: this

    ! Kill and deallocate testTransportNuclearData
    call this % xsData % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test correct behaviour of the filter
  !!
@Test
  subroutine testGettingResponse(this)
    class(test_macroResponse), intent(inout) :: this
    type(particle)                           :: p
    real(defReal), parameter :: TOL = 1.0E-9

    p % type = P_NEUTRON
    p % isMG = .false.

    ! Test response values
    @assertEqual(6.0_defReal, this % response_total % get(p, this % xsData), TOL)
    @assertEqual(2.0_defReal, this % response_capture % get(p, this % xsData), TOL)
    @assertEqual(1.0_defReal, this % response_fission % get(p, this % xsData), TOL)
    @assertEqual(1.5_defReal, this % response_nuFission % get(p, this % xsData), TOL)
    @assertEqual(3.0_defReal, this % response_absorbtion % get(p, this % xsData), TOL)
    @assertEqual(105.0_defReal, this % response_MT % get(p, this % xsData), TOL)

    p % isMG = .true.
    @assertEqual(ZERO, this % response_MT % get(p, this % xsData), TOL)

  end subroutine testGettingResponse

end module macroResponse_test
