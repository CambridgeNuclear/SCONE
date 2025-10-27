module microResponse_test

  use numPrecision
  use endfConstants
  use microResponse_class,            only : microResponse
  use particle_class,                 only : particle, P_NEUTRON
  use dictionary_class,               only : dictionary
  use testNeutronDatabase_class,      only : testNeutronDatabase
  use materialMenu_mod,               only : init
  use funit

  implicit none

@testCase
  type, extends(TestCase) :: test_microResponse
    private
    type(microResponse)        :: response_total
    type(microResponse)        :: response_eScatter
    type(microResponse)        :: response_capture
    type(microResponse)        :: response_fission
    type(microResponse)        :: response_absorbtion
    type(microResponse)        :: response_MT
    type(testNeutronDatabase)  :: xsData
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_microResponse


contains

  !!
  !! Sets up test_microResponse object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_microResponse), intent(inout) :: this
    type(dictionary)                         :: tempDict, dictMat1, dictMat2, dictMat3

    ! Allocate and initialise test nuclearData

    ! Cross-sections:         Total         eScattering  IeScatter  Capture     Fission       nuFission
    call this % xsData % build(6.0_defReal, 3.0_defReal, ZERO,     2.0_defReal, 1.0_defReal, 1.5_defReal)

    ! Set dictionaries to initialise material
    call dictMat1 % init(1)
    call dictMat2 % init(2)
    call dictMat3 % init(1)

    call dictMat3 % store('54135.03', 2.0_defReal)

    call dictMat2 % store('temp', 300.0_defReal)
    call dictMat2 % store('composition', dictMat3)

    call dictMat1 % store('Xenon', dictMat2)

    ! Initialise material
    call init(dictMat1)

    ! Set up responses
    ! Total
    call tempDict % init(3)
    call tempDict % store('type', 'microResponse')
    call tempDict % store('MT', N_TOTAL)
    call tempDict % store('material', 'Xenon')
    call this % response_total % init(tempDict)
    call tempDict % kill()

    ! Capture
    call tempDict % init(3)
    call tempDict % store('type', 'microResponse')
    call tempDict % store('MT', N_DISAP)
    call tempDict % store('material', 'Xenon')
    call this % response_capture % init(tempDict)
    call tempDict % kill()

    ! Fission
    call tempDict % init(3)
    call tempDict % store('type', 'microResponse')
    call tempDict % store('MT', N_FISSION)
    call tempDict % store('material', 'Xenon')
    call this % response_fission % init(tempDict)
    call tempDict % kill()

    ! nuFission
    call tempDict % init(3)
    call tempDict % store('type', 'microResponse')
    call tempDict % store('MT', N_N_ELASTIC)
    call tempDict % store('material', 'Xenon')
    call this % response_eScatter % init(tempDict)
    call tempDict % kill()

    ! Absorbtion
    call tempDict % init(3)
    call tempDict % store('type', 'microResponse')
    call tempDict % store('MT', N_ABSORPTION)
    call tempDict % store('material', 'Xenon')
    call this % response_absorbtion % init(tempDict)
    call tempDict % kill()

    ! MT number
    call tempDict % init(2)
    call tempDict % store('type','microResponse')
    call tempDict % store('MT', 16)
    call tempDict % store('material', 'Xenon')
    call this % response_MT % init(tempDict)
    call tempDict % kill()

  end subroutine setUp

  !!
  !! Kills test_microResponse object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_microResponse), intent(inout) :: this

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
    class(test_microResponse), intent(inout) :: this
    type(particle)                           :: p
    real(defReal), parameter :: TOL = 1.0E-9

    p % type = P_NEUTRON
    p % isMG = .false.

    ! Test response values
    @assertEqual(3.0_defReal, this % response_total % get(p, this % xsData), TOL)
    @assertEqual(1.0_defReal, this % response_capture % get(p, this % xsData), TOL)
    @assertEqual(0.5_defReal, this % response_fission % get(p, this % xsData), TOL)
    @assertEqual(1.5_defReal, this % response_eScatter % get(p, this % xsData), TOL)
    @assertEqual(1.5_defReal, this % response_absorbtion % get(p, this % xsData), TOL)
    @assertEqual(8.0_defReal, this % response_MT % get(p, this % xsData), TOL)

    p % isMG = .true.
    @assertEqual(ZERO, this % response_MT % get(p, this % xsData), TOL)

  end subroutine testGettingResponse

end module microResponse_test
