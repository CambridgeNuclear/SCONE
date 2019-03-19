!!
!! Timer module
!!
!! Is effectively a singleton to store information about the runtime
!! It is good to think about this module as a global class.
!!
!! Interface:
!!   stopWatch -> class that can be used to time elapsed time. See its docs for details
!!   secToChar -> conversion function. Takes number of seconds as a real. Returns character
!!                of length nameLen with left-adjusted time in hhh:mm:ss format
!!
!!
!! Private members:
!!
!!
!!
module timer_mod

  use numPrecision

  implicit none
  private

  !! Public interface
  public :: secToChar

  !!
  !! This derived type allows to measure time in a program
  !!
  !! It is using system_clock to have a best possible resolution for very short
  !! procedures. Given the size of 64-bit int and only nano-second resolution of the
  !! system clock problems with an overflow are unlikley to happen
  !!
  !! Members:
  !!  start_c   -> count on start
  !!  elapsed_c -> total accumulated time in clock counts
  !!
  !! Interface:
  !!   start  -> start stopWatch
  !!   stop   -> end counting and add result to elapsed time
  !!   reset  -> set accumulated elapsed time to 0
  !!   toSec  -> returns real number with the estimate of seconds accumulated
  !!   toChar -> print elapsed time to human readable format hhh:mm:ss format
  !!
  type,public :: stopWatch
    private
    integer(longInt) :: start_c   = 0
    integer(longInt) :: elapsed_c = 0
  contains
    procedure :: start  => start_stopWatch
    procedure :: stop   => stop_stopWatch
    procedure :: reset  => reset_stopWatch
    procedure :: toSec  => toSec_stopWatch
  end type stopWatch

contains

  !!
  !! Return a nameLen string with elapsed time in hhh:mm:ss format
  !!
  elemental function secToChar(sec) result(time)
    real(defReal),intent(in)     :: sec
    real(defReal)                :: elapsed
    character(nameLen)           :: time
    integer(shortInt)            :: hour, minute, second

    ! Copy number of seconds into work variable
    elapsed = sec

    ! Calculate number of hours
    hour = floor(elapsed / 3600.0_defReal)

    ! Calculate number of minutes
    elapsed = elapsed - hour * 3600.0_defReal
    minute  = floor(elapsed / 60.0_defReal)

    ! Calculate number of seconds
    elapsed = elapsed - minute * 60.0_defReal
    second = floor(elapsed)

    write(time,'(I3,A1,I2.2,A1,I2.2 )') hour, ':', minute, ':', second

  end function secToChar


!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! StopWatch Procedures
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Start stopWatch
  !! Sets new value into start memeber
  !!
  subroutine start_stopWatch(self)
    class(stopWatch), intent(inout) :: self

    call system_clock(count = self % start_c)

  end subroutine start_stopWatch

  !!
  !! Stops stopWatch
  !! Adds counts to the elapsed time
  !! Replaces start with current count
  !!
  subroutine stop_stopWatch(self)
    class(stopWatch), intent(inout) :: self
    integer(longInt)                :: ticks

    ticks = self % start_c
    call system_clock(count = self % start_c)
    self % elapsed_c = self % elapsed_c + self % start_c - ticks

  end subroutine stop_stopWatch

  !!
  !! Reset elapsed time to 0
  !!
  subroutine reset_stopWatch(self)
    class(StopWatch), intent(inout) :: self

    self % elapsed_c = 0

  end subroutine reset_stopWatch

  !!
  !! Return elapsed time in seconds
  !!
  function toSec_stopWatch(self) result(sec)
    class(stopWatch), intent(in) :: self
    real(defReal)                :: sec
    integer(longInt)             :: rate

    call system_clock(count_rate = rate)
    sec = real(self % elapsed_c, defReal) / rate

  end function toSec_stopWatch

end module timer_mod
