module imgBmp_iTest

  use numPrecision
  use imgBmp_func, only : intToByte, imgBmp
  use pFUnit_mod

  implicit none


contains

  !!
  !! Test little endiness byte printng
  !!
@Test
  subroutine test_byte_print()
    character(1) :: ref1
    character(2) :: ref2
    character(3) :: ref3
    character(4) :: ref4

    ! 1 Byte value
    ref1 = transfer(z'86', ref1)   ! 134

    ! 2 Byte value -> Construct byte by byte Necessary to preserve endianess
    ! 259
    ref2(1:1) = transfer(z'03', ref1)
    ref2(2:2) = transfer(z'01', ref1)

    ! 3 Byte value -
    ! 776655
    ref3(1:1) = transfer(z'CF', ref1)
    ref3(2:2) = transfer(z'D9', ref1)
    ref3(3:3) = transfer(z'0B', ref1)

    ! 4 Byte Value
    ! 133316666
    ref4(1:1) = transfer(z'3A', ref1)
    ref4(2:2) = transfer(z'40', ref1)
    ref4(3:3) = transfer(z'F2', ref1)
    ref4(4:4) = transfer(z'07', ref1)

    ! Verify
    @assertEqual(ref1, intToByte(134, 1))
    @assertEqual(ref2, intToByte(259, 2))
    @assertEqual(ref3, intToByte(776655, 3))
    @assertEqual(ref4, intToByte(133316666, 4))

  end subroutine test_byte_print

  !!
  !! Test creation of a BMP image by comparing the output to the reference BMP file
  !!
@Test
  subroutine test_image()
    integer(shortInt), dimension(5,5) :: img
    integer(shortInt)                 :: black, green, red, yellow, white, blue
    integer(shortInt)                 :: file, i
    character(1)                      :: ref
    character(:), allocatable         :: image
    character(*), parameter           :: path = './IntegrationTestFiles/sample.bmp'

    ! Write colors
    black  = int(z'000000', shortInt)
    white  = int(z'FFFFFF', shortInt)
    yellow = int(z'FFFF00', shortInt)
    green  = int(z'00FF00', shortInt)
    red    = int(z'0000FF', shortInt)
    blue   = int(z'FF0000', shortInt)

    img = black
    img(1,1) = white
    img(5,1) = yellow
    img(3,3) = green
    img(1,5) = red
    img(5,5) = blue

    ! Create Image as a character string
    image = imgBmp(img)

    ! Compare byte by byte with the reference
    open(newunit=file, file=path, access='stream', status='old', action='read')

    do i = 1, len(image)
      read(file) ref
      @assertEqual(ref, image(i:i))
    end do

    close(file)

  end subroutine test_image

end module imgBmp_iTest
