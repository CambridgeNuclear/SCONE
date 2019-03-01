module testTransportNuclearData_test

  use numPrecision
  use testTransportNuclearData_class, only : testTransportNuclearData
  use xsMacroSet_class,               only : xsMacroSet_ptr
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

    call this % nucData % build(2.0_defReal,              &
                                scatterXS   = 3.0_defReal,&
                                captureXS   = 4.0_defReal,&
                                fissionXS   = 5.0_defReal,&
                                nuFissionXS = 6.0_defReal)

  end subroutine setUp

  !!
  !! Kills test_testTransportNuclearData object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_testTransportNuclearData), intent(inout) :: this

    call this % nucData % kill()

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

  !!
  !! Test obtaining xsMacroSet
  !!
@Test
  subroutine testMacroXsSet(this)
    class(test_testTransportNuclearData), intent(inout) :: this
    type(particle)                                      :: p
    type(testTransportNuclearData)                      :: nucDat
    type(xsMacroSet_ptr)                                :: XSs
    real(defReal),parameter                             :: TOL = 1.0E-9

    ! Get xsMacroSet ptr
    call this % nucData % getMatMacroXS(XSs, p, 1)

    @assertEqual(2.0_defReal,XSs % totalXS(), TOL)
    @assertEqual(3.0_defReal,XSs % scatterXS(), TOL)
    @assertEqual(4.0_defReal,XSs % captureXS(), TOL)
    @assertEqual(5.0_defReal,XSs % fissionXS(), TOL)
    @assertEqual(6.0_defReal,XSs % nuFissionXS(), TOL)

    ! Check with diffrent initialisation
    call nucDat % build(ONE)

    call nucDat % getMatMacroXS(XSs, p, 1)
    @assertEqual(ONE, XSs % totalXS(), TOL)
    @assertEqual(ONE, XSs % scatterXS(), TOL)
    @assertEqual(ONE, XSs % captureXS(), TOL)
    @assertEqual(ONE, XSs % fissionXS(), TOL)
    @assertEqual(ONE, XSs % nuFissionXS(), TOL)

    call nucDat % kill()

  end subroutine testMacroXsSet

  !!
  !! Test isFissile
  !!
@Test
  subroutine testIsFissile(this)
    class(test_testTransportNuclearData), intent(inout) :: this
    type(testTransportNuclearData)                      :: nucDat

    @assertTrue(this % nucData % isFissileMat(1))
    call nucDat % build(ONE, fissionXS = ZERO)
    @assertFalse(nucDat % isFissileMat(1))

    call nucDat % kill()
  end subroutine testIsFissile

end module testTransportNuclearData_test
