!!
!! Functions that allow printing elementary bitmap images
!!
!! Only 24-bit Color Bitmaps (without any compression or color palette data) can be created.
!!
!! It is not ideal due to large size of the BMP files. However, the purpose of theese functions
!! is to provide some image creataion capability without external dependencies.
!!
!! Implementation should work even when the Fortran compiler prints data in big-endian format.
!! (I guess it depends on the compiler+system combination. I do wander what gfortran would
!!  do for the bi-endian CPU)
!!
!! Has been based on:
!!   https://itnext.io/bits-to-bitmaps-a-simple-walkthrough-of-bmp-image-format-765dc6857393
!!
!! Interface:
!!   imgBmp        -> Generate bitmap image as a string of characters (bytes)
!!   imgBmp_toFile -> Generate bitmap image as a file
!!
module imgBmp_func

  use numPrecision
  use genericProcedures, only : fatalError, numToChar

  implicit none
  private

  ! Local parameters
  ! Get kind of char to ensure it is exactly a Byte
  integer(shortInt), parameter :: byte = selected_char_kind('ascii')

  ! Public Interface
  public :: intToByte ! Public for testing
  public :: imgBmp
  public :: imgBmp_toFile


contains

  !!
  !! Print BMP image directly to file
  !!
  !! See `imgBmp` for details on the BMP requirements.
  !!
  !! Args:
  !!   map [in] -> Integer matrix with pixel data. See `imgBmp` for more details
  !!   file [in] -> Character with path to the output file
  !!
  !! Errors:
  !!   If path in `file` points to an existing file overrride it.
  !!
  subroutine imgBmp_toFile(map, file)
    integer(shortInt), dimension(:,:), intent(in) :: map
    character(*), intent(in)                      :: file
    integer(shortInt)                             :: unit

    open(newunit=unit, file=file, status='replace', access='stream', action='write')

    write(unit) imgBmp(map)

    close(unit)

  end subroutine imgBmp_toFile

  !!
  !! Print BMP image to string of bytes
  !!
  !! Creates a 24-bit color bitmap from the integer array.
  !! Does not allow for colormaps, compression or any advanced features.
  !! Map must contain integers representing 24-bit colors. Columns of map
  !! represents rows of the bitmap. Thus map(3,4) is the 3rd pixel right and 4th above
  !! the lower left corner.
  !!
  !! Args:
  !!   map [in] -> Integer NxM array with color for each pixel.
  !!
  !! Result:
  !!   String of bytes that is a binary representation of the BMP image.
  !!
  !! Errors:
  !!   fatalError if map contains invalid integers (-ve or larger then 24-bits)
  !!
  function imgBmp(map) result(str)
    integer(shortInt), dimension(:,:), intent(in) :: map
    character(:), allocatable                     :: str
    integer(shortInt)                             :: offset, L, padBits, NX, NY
    character(100), parameter :: Here = 'imgBmp (imgBmp_func.f90)'

    ! Check map
    if (any(map < 0 .or. map >= 2**24 )) then
      call fatalError(Here, 'Image integer map contains values that are -ve or larger &
                            &then 24-bit integer.')
    end if

    ! Get sizes for conveniance
    NX = size(map, 1)
    NY = size(map, 2)
    padBits = modulo(4 - modulo(3*NX, 4), 4) ! A bit of a one-liner for now

    ! Calculate total size
    L = 14 + 40 + (3*NX + padBits) * NY
    offset = 14 + 40

    ! Create image bit-string
    str = bmpHead(L, offset) // bmpImgInfo(NX, NY) // bmpPixelData(map, padBits)

  end function imgBmp

  !!
  !! Print unsigned integer to bytestring in little endian layout
  !!
  !! Should also work on Big
  !!
  !! Args:
  !!   num [in] -> Positive integer (given as signed)
  !!   N [in]   -> Number of bytes for the integer [1-4]
  !!
  !! Result:
  !!   N bytes representing the integer in little endian layout.
  !!
  !! Errors:
  !!   fatalError if num requires more then N bytes.
  !!   fatalError if the input number is -ve
  !!
  function intToByte(num, N) result(str)
    integer(shortInt), intent(in) :: num
    integer(shortInt), intent(in) :: N
    character(N, byte)            :: str
    integer(shortInt)             :: i, temp, max
    integer(shortInt), parameter  :: mask8 = 255
    character(100), parameter :: Here = 'intToByte (imgBmp_func.f90)'

    ! Generate Maximum Integer that can be represented
    if ( N <= 3) then
      max = ibset(0, N*8)
    else
      max = huge(num)
    end if

    ! Check for errors
    if (num < 0) then
      call fatalError(Here, 'Was given -ve number: '//numToChar(num))

    else if (num >= max) then
      call fatalError(Here, 'Integer: '//numToChar(num)//' does not fit into '//&
                            numToChar(N)//' bytes.')
    end if

    ! Convert
    ! By elier check we can assert that sign-bit (most significant bit) is 0
    temp = num
    do i = 1, 4
      str(i:i) = char(iand(temp, mask8), byte)
      temp = temp / 256
    end do

  end function intToByte

  !!
  !! Print BMP Header
  !!
  !! Args:
  !!   N [in] -> Total size of the BMP file in bytes
  !!   offset [in] -> Number of bytes between beginning of the file and start of the image data
  !!     (total size of the header, info and color palette blocks)
  !!
  !! Result:
  !!   String of 14 bytes with the header
  !!
  function bmpHead(N, offset) result(str)
    integer(shortInt), intent(in) :: N
    integer(shortInt), intent(in) :: offset
    character(14, byte)           :: str
    character(1, byte)            :: mold

    ! Print file type
    str(1:1) = transfer(z'42', mold)
    str(2:2) = transfer(z'4D', mold)

    ! Print size
    str(3:6) = intToByte(N, 4)

    ! Few empty bytes
    str(7:8) = intToByte(0, 2)
    str(9:10) = intToByte(0, 2)

    !
    str(11:14) = intToByte(offset, 4)

  end function bmpHead

  !!
  !! Print Bitmap Information Header (DIB)
  !!
  !! Prints the 40-byte DIB block for the bitmap image that uses 24-bit colors.
  !!
  !! Args:
  !!   width [in] -> Width (x-axis) of the image in pixels
  !!   height [in] -> Hight (y-axis) of the image in pixels
  !!
  !! Result:
  !!   40 byte string with DIB information. Bare bones version only with
  !!   no compression nor color palette. Does not set the resolution (pixels/m) data
  !!
  function bmpImgInfo(width, height) result(str)
    integer(shortInt), intent(in) :: width
    integer(shortInt), intent(in) :: height
    character(40, byte)           :: str

    ! Print header size
    str(1:4) = intToByte(40, 4)

    ! Print width & height
    str(5:8)  = intToByte(width, 4)
    str(9:12) = intToByte(height, 4)

    ! Planes of target device
    str(13:14) = intToByte(1, 2)

    ! Bits per pixel
    str(15:16) = intToByte(24, 2)

    ! Compresion info
    str(17:20) = intToByte(0, 4)
    str(21:24) = intToByte(0, 4)

    ! Pixels per metere
    str(25:28) = intToByte(0, 4)
    str(29:32) = intToByte(0, 4)

    ! Color Info -> No color map
    str(33:36) = intToByte(0, 4)
    str(37:40) = intToByte(0, 4)

  end function bmpImgInfo

  !!
  !! Print BMP Pixel Data
  !!
  !! Print BMP image pixel data from NxM array of 24-bit colors.
  !!
  !! Args:
  !!   map [in] -> Integer NxM array with color for each pixel. Columns (1st index) of the
  !!     matrix represent rows. map(3,4) is the 3rd pixel right and 4th above the lower left corner.
  !!   padbyte [in] -> Required number of padding bits (0-3) so the row data size is multiple of 4
  !!     bytes.
  !!
  !! Result:
  !!   String with the image pixel data.
  !!
  !! Errors:
  !!   fatalError if padding is outside the valid range
  !!
  function bmpPixelData(map, padbyte) result(str)
    integer(shortInt), dimension(:,:), intent(in) :: map
    integer(shortInt), intent(in)                 :: padbyte
    character(:), allocatable                     :: str
    integer(shortInt)                             :: NX, NY, N_total, pos, i, j
    character(100), parameter :: Here = 'bmpPixelData (imgBmp_func.f90)'

    ! Check padding
    if ( padbyte < 0 .or. 3 < padbyte) then
      call fatalError(Here, 'Invalid number of padding bytes. &
                            &Must be 1-3, is: '//numToChar(padbyte))
    end if

    ! Calculate total size
    NX = size(map, 1)
    NY = size(map, 2)
    N_total = (NX * 3 + padbyte) * NY

    allocate(character(N_total, byte) :: str)

    ! Print colors
    pos = 1
    do j = 1, NY
      do i = 1, NX
        ! Print color
        str(pos:pos+2) = intToByte(map(i,j), 3)
        pos = pos + 3
      end do
      ! Print padding
      str(pos:pos+padbyte-1) = intToByte(0, padbyte)
      pos = pos + padbyte

    end do

  end function bmpPixelData

end module imgBmp_func
