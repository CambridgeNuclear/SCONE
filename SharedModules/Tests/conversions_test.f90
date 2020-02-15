module conversions_test
  use numPrecision
  use genericProcedures, only : charToInt
  use pFUnit_mod

  implicit none

contains

  !!
  !! Test char to int conversion
  !!
@Test
  subroutine testCharToInt()
    logical(defBool)  :: flag
    integer(shortInt) :: i
    ! Easy cases
    @assertEqual(2, charToInt('   2 '))
    @assertEqual(-1, charToInt('  -1'))

    ! Edge cases
    @assertEqual(7,charToInt('7                   A'))


    ! Cases that should fail
    i = charToInt('7.0', error = flag)
    @assertTrue(flag)

    i = charToInt('7E-7', error = flag)
    @assertTrue(flag)

    i = charToInt('NaN', error = flag)
    @assertTrue(flag)

    i = charToInt('7 is Not A Number For Fortran Rex', error = flag)
    @assertTrue(flag)

    i = charToInt('253.', error = flag)
    @assertTrue(flag)

  end subroutine testCharToInt

    
end module conversions_test
