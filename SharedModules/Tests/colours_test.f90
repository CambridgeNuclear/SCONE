module colours_test
  use numPrecision
  use colours_func, only : rgb24bit
  use funit

  implicit none

contains


@Test
  subroutine testColourConversions()

    ! Test by comparison with some hex values
    @assertEqual(int(z"123456"), rgb24bit(18, 52, 86))

    @assertEqual(int(z"000000"), rgb24bit(0, 0, 0))
    @assertEqual(int(z"ffffff"), rgb24bit(255, 255, 255))

    @assertEqual(int(z"ff0000"), rgb24bit(255, 0, 0))
    @assertEqual(int(z"00ff00"), rgb24bit(0, 255, 0))
    @assertEqual(int(z"0000ff"), rgb24bit(0, 0, 255))

  end subroutine testColourConversions

end module colours_test
