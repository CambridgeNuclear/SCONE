!!
!! Contains function to deal with colours
!!
module colours_func
  use numPrecision
  use genericProcedures, only: numToChar
  use errors_mod,        only: fatalError
  implicit none

  public :: rgb24bit

contains


  !!
  !! Return integer with 24bit colour value from r,g,b components
  !!
  !! Args:
  !!   r [in] -> red component in [0-255] range
  !!   g [in] -> green component in [0-255] range
  !!   b [in] -> blue component in [0-255] range
  !!
  function rgb24bit(r, g, b) result(col)
    integer(shortInt), intent(in) :: r
    integer(shortInt), intent(in) :: g
    integer(shortInt), intent(in) :: b
    integer(shortInt)             :: col
    character(100), parameter :: Here = "rgb24bit (colours_func.f90)"

    ! Check that the values are in the correct range
    if (r < 0 .or. r > 255) then
      call fatalError(Here, "Red value is out of [0-255] range: " // numToChar(r))
    end if
    if (g < 0 .or. g > 255) then
      call fatalError(Here, "Green value is out of [0-255] range: " // numToChar(g))
    end if
    if (b < 0 .or. b > 255) then
      call fatalError(Here, "Blue value is out of [0-255] range: " // numToChar(b))
    end if

    ! Compute the 24 bit colour
    ! Note inversion, blue is the least significant bits
    col = b + 256 * g + 256 * 256 * r

  end function rgb24bit


end module colours_func
