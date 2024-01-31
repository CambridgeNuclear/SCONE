module delayedStream_class

  use numPrecision
  use errors_mod, only : fatalError

  implicit none
  private

  integer(shortInt), parameter :: MAX_DELAY = 1024
  integer(shortInt), parameter :: MIN_DELAY = 16

  !!
  !! Writes characters to a file with a delay to allow for backtracking.
  !!
  !! It is required to make writing of output format in stream-like fashion
  !! easier to implement. For example, in JSON, it allows to backtrack and
  !! remove trailing comma when closing a JSON dictionary.
  !!
  !! NOTE: We need to be careful to check if the file is open before writing
  !!       to it. Otherwise, we will get an error.
  !!
  !! Private Member:
  !!   unit      -> Unit to which the output file is connected
  !!   bufferPos -> Position in the buffer
  !!   buffer    -> Buffer of characters
  !!
  !! Interface:
  !!   write -> Write some characters to the buffer/file
  !!   cut   -> Remove last N characters from the buffer
  !!   peek  -> Inspect the Nth last character in the buffer
  !!   close -> Write all characters remaining in the buffer to the file
  !!   setUnit -> Set the unit to which the output file is connected
  !!
  type, public :: delayedStream
    private
    integer(shortInt) :: unit = -8
    integer(shortInt) :: bufferPos = 0
    character(MAX_DELAY) :: buffer
  contains
    procedure :: write
    procedure :: cut
    procedure :: peek
    procedure :: close
    procedure :: setUnit

    procedure, private :: flush
  end type delayedStream

contains

  !!
  !! Write to a delayed stream
  !!
  !! Args:
  !!   text [in] -> text to write
  !!
  subroutine write(self, text)
    class(delayedStream), intent(inout) :: self
    character(*), intent(in) :: text
    integer(shortInt) :: i

    ! Write to buffer
    do i = 1, len(text)
      self % bufferPos = self % bufferPos + 1
      self % buffer(self % bufferPos:self % bufferPos) = text(i:i)
      if (self % bufferPos == MAX_DELAY) then
        call flush(self)
      end if
    end do

  end subroutine write

  !!
  !! Remove last characters from a delayed stream
  !!
  !! Args:
  !!  n [in] -> number of characters to remove. Must be less than the MIN_DELAY.
  !!            If n is larger than the number of characters in buffer, everything is removed.
  !!
  subroutine cut(self, n)
    class(delayedStream), intent(inout) :: self
    integer(shortInt), intent(in) :: n
    character(100), parameter :: HERE = "cut (delayedStream_class.f90)"

    if (n > MIN_DELAY) then
      call fatalError(HERE, "n must be less than MIN_DELAY")
    end if

    self % bufferPos = max(0, self % bufferPos - n)

  end subroutine cut

  !!
  !! Inspect the nth last character in the buffer
  !!
  !! Args:
  !!  n [in] -> The number of characters to peek back.
  !!
  function peek(self, n) result(c)
    class(delayedStream), intent(inout) :: self
    integer(shortInt), intent(in) :: n
    character(n) :: c
    character(100), parameter :: HERE = "peek (delayedStream_class.f90)"

    if (n > self % bufferPos) then
      call fatalError(HERE, "n is larger than the buffer")
    else
      c = self % buffer(self % bufferPos - n + 1: self % bufferPos - n + 1)
    end if

  end function peek

  !!
  !! Flush contents of the buffer
  !!
  subroutine flush(self)
    class(delayedStream), intent(inout) :: self
    integer(shortInt) :: to_file
    logical(defBool) :: isOpen

    ! Do not write anything if the file is not opened
    inquire(unit=self % unit, opened=isOpen)
    if (.not. isOpen) then
      return
    end if

    ! Do nothing if the buffer is below minimum length
    if (self % bufferPos <= MIN_DELAY) then
      return
    end if
    to_file = self % bufferPos - MIN_DELAY

    ! Write to file only if it is open
    if (isOpen) then
      write(self % unit, "(A)", advance="no") self % buffer(1 : to_file)
    end if

    self % buffer(1:MIN_DELAY) = self % buffer(to_file + 1 : self % bufferPos)

    self % bufferPos = MIN_DELAY

  end subroutine flush

  !!
  !! Close the delayed stream
  !!
  subroutine close(self)
    class(delayedStream), intent(inout) :: self
    logical(defBool) :: isOpen

    ! Do not write anything if the file is not opened
    inquire(unit=self % unit, opened=isOpen)
    if (.not. isOpen) then
      return
    end if

    ! Does nothing if empty
    if (self % bufferPos == 0) then
      return
    end if

    ! Write the remaining buffer
    ! Write to file only if it is open
    if (isOpen) then
      write(self % unit, "(A)", advance="no") self % buffer(1:self % bufferPos)
    end if
    self % bufferPos = 0

  end subroutine close

  !!
  !! Set the unit to which the output file is connected
  !!
  !! Args:
  !!  unit [in] -> Unit to which the output file is connected. May be opened or not.
  !!               Both are fine. If it is not opened, nothing will be written to it.
  subroutine setUnit(self, unit)
    class(delayedStream), intent(inout) :: self
    integer(shortInt), intent(in) :: unit

    self % unit = unit

  end subroutine setUnit

end module delayedStream_class
