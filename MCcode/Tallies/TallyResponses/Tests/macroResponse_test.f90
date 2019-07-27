module macroResponse_test

  use numPrecision
  use endfConstants
  use macroResponse_class,            only : macroResponse
  use particle_class,                 only : particle
  use dictionary_class,               only : dictionary
  use testTransportNuclearData_class, only : testTransportNuclearData
  use pFUnit_mod

  implicit none

@testCase
  type, extends(TestCase) :: test_macroResponse
    private
    type(macroResponse)                     :: response_total
    type(macroResponse)                     :: response_capture
    type(macroResponse)                     :: response_fission
    type(macroResponse)                     :: response_nuFission
    type(macroResponse)                     :: response_absorbtion
    type(testTransportNuclearData), pointer :: xsData => null()
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
    allocate(this % xsData)

    ! Cross-sections:         Total        Scatering    Capture     Fission       nuFission
    call this % xsData % build(6.0_defReal, 3.0_defReal, 2.0_defReal, 1.0_defReal, 1.5_defReal)

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

  end subroutine setUp

  !!
  !! Kills test_macroResponse object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_macroResponse), intent(inout) :: this

    ! Kill and deallocate testTransportNuclearData
    call this % xsData % kill()
    deallocate(this % xsData)

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

    ! Prepare particle
    p % xsData => this % xsData

    ! Test response values
    @assertEqual(6.0_defReal, this % response_total % get(p), TOL)
    @assertEqual(2.0_defReal, this % response_capture % get(p), TOL)
    @assertEqual(1.0_defReal, this % response_fission % get(p), TOL)
    @assertEqual(1.5_defReal, this % response_nuFission % get(p), TOL)
    @assertEqual(3.0_defReal, this % response_absorbtion % get(p), TOL)

  end subroutine testGettingResponse

end module macroResponse_test
