module testTransportNuclearData_test

  use numPrecision
  use testTransportNuclearData_class, only : testTransportNuclearData
  use particle_class,                 only : particle
  use pFUnit_mod

  implicit none

@testCase
  type, extends(TestCase) :: test_testTransportNuclearData
    private
    type(testTransportNuclearData) :: nucData
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_testTransportNuclearData

contains

  !!
  !! Sets up test_testTransportNuclearData object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_testTransportNuclearData), intent(inout) :: this

    call this % nucData % build(2.0_defReal)

  end subroutine setUp

  !!
  !! Kills test_testTransportNuclearData object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_testTransportNuclearData), intent(inout) :: this

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test getting totalXS
  !!
@Test
  subroutine testGetXSs(this)
    class(test_testTransportNuclearData), intent(inout) :: this
    type(particle)                                      :: p
    real(defReal)                                       :: xsVal
    real(defReal)                                       :: tol = 1.0E-9

    ! Get trans XS
    xsVal = this % nucData % getTransXS(p, 1)
    @assertEqual(2.0_defReal, xsVal, tol)

    ! Get majorant XS
    xsVal = this % nucData % getMajorantXS(p)
    @assertEqual(2.0_defReal, xsVal, tol)

    ! Get total XS
    xsVal = this % nucData % getTotalMatXS(p, 1)
    @assertEqual(2.0_defReal, xsVal, tol)

  end subroutine testGetXSs

end module testTransportNuclearData_test
