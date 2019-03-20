!!
!! Timer module
!!
!! Is effectively a singleton to store information about the runtime
!! It is good to think about this module as a global class.
!!
!! TODO: Make provisions for OpenMP safety.
!!
!!
!! Interface:
!!   stopWatch -> class that can be used to time elapsed time. See its docs for details
!!   secToChar -> conversion function. Takes number of seconds as a real. Returns character
!!                of length nameLen with left-adjusted time in hhh:mm:ss format
!!   registerTimer -> create new timer definition.
!!   timerStart    -> Start a defined timer
!!   timerStop     -> Stop a defined timer
!!   timerReset    -> Reset a defined timer
!!   timerTime     -> Get total elapsed time in Seconds
!!
!! Private members:
!!   timerNames -> array of names of different timers
!!   timers     -> array of stopwatches
!!   idx        -> index of the next avalible bin
!!
!!
module timer_mod

  use numPrecision
  use genericProcedures, only : fatalError

  implicit none
  private

  !! Public interface
  public :: secToChar
  public :: registerTimer
  public :: timerStart
  public :: timerStop
  public :: timerReset
  public :: timerTime

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

  !! Module members
  character(nameLen),dimension(:),allocatable :: timerNames
  type(stopWatch),dimension(:),allocatable    :: timers
  integer(shortInt)                           :: idx = 0

  !! Module parameters
  real(defReal),parameter     :: GROWTH_RATIO = 1.6_defReal
  integer(shortInt),parameter :: MIN_SIZE = 5

contains

  !!
  !! Ensure that there is enough memory for new timer
  !! Grow and reallocate space if needed
  !!
  !! Args:
  !!   None
  !!
  !! Errors:
  !!   Returns fatalError if by evil magic allocation status of "timers" and "timerNames" is
  !!   diffrent. This should never happen.
  !!
  !! Notes:
  !!   Uses "idx" to determine if reallocation is needed.
  !!   Reallocates if idx > current size
  !!   For unallocated "timers" and "timerNames" current size is set to 0!
  !!   "idx" needs to be increamented BEFORE checking for memorty growth.
  !!
  subroutine growIfNeeded()
    logical(defBool)                            :: nameAlloc, timersAlloc
    integer(shortInt)                           :: newSize, S
    character(nameLen),dimension(:),allocatable :: namesTemp
    type(stopWatch),dimension(:),allocatable    :: timerTemp
    character(100),parameter :: Here ='growIfNeeded (timer_mod.f90)'

    ! Check allocation status
    nameAlloc   = allocated(timerNames)
    timersAlloc = allocated(timers)

    ! Set current size for unallocated case
    S = 0

    ! Calculate the required size
    if(nameAlloc .and. timersAlloc) then
      ! If allocated load real size
      S = size(timers)
      newSize = S

      ! Check if reallocation is needed and calculate new size
      if (idx >= size(timers)) newSize = int(newSize * GROWTH_RATIO)

    else if( nameAlloc .neqv. timersAlloc) then ! Strange error has happend!
      call fatalError(Here,'Timers and timers names array have diffrent allocation status. WTF?')
      newSize = 0

    else ! Unallocated case. Allocate to minimum size.
      newSize = MIN_SIZE

    end if

    ! Reallocate if needed
    if( newSize > S) then
      allocate(namesTemp(newSize))
      allocate(timerTemp(newSize))

      ! Copy only if already allocated. Avoid SEG errors
      if(nameAlloc)   namesTemp(1:idx) = timerNames
      if(timersAlloc) timerTemp(1:idx) = timers

      ! Move allocation
      call move_alloc(namesTemp, timerNames)
      call move_alloc(timerTemp, timers)
    end if
  end subroutine growIfNeeded

  !!
  !! Return a nameLen string with elapsed time in hhh:mm:ss format
  !!
  !! Args:
  !!  sec [in] -> real(defReal) number with number of second
  !!
  !! Returns:
  !!   nameLen long character with left-adjusted time in hhh:mm:ss format
  !!
  !! Errors:
  !!   None.
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

  !!
  !! Register timer Function
  !!
  !! Creates new bin for the timer.
  !!
  !! Args:
  !!   name [in] -> new timer name. Character of maximum length of nameLen
  !!
  !! Returns:
  !!   Integer with bin index
  !!
  !! Errors:
  !!   None
  !!
  function registerTimer(name) result(timer_idx)
    character(*), intent(in) :: name
    integer(shortInt)        :: timer_idx

    ! Increment count
    idx = idx + 1
    timer_idx = idx

    ! Grow storage space if needed
    call growIfNeeded()

    ! Load name for the timer
    timerNames(idx) = name

  end function registerTimer

  !!
  !! timerStart Subroutine
  !!
  !! Start timer under index
  !!
  !! Args:
  !!   binIdx [in] - Index of the timer, recived from registerTimer function
  !!
  !! Errors:
  !!   If binIdx is out-of-bounds, behaviour is udefined. Most likely SEG error
  !!
  subroutine timerStart(binIdx)
    integer(shortInt), intent(in) :: binIdx

    call timers(binIdx) % start()

  end subroutine timerStart

  !!
  !! timerStop Subroutine
  !!
  !! Stops timer under index. Add counts to total elapsed time
  !!
  !! Args:
  !!   binIdx [in] - Index of the timer, recived from registerTimer function
  !!
  !! Errors:
  !!   If binIdx is out-of-bounds, behaviour is udefined. Most likely SEG error
  !!
  subroutine timerStop(binIdx)
    integer(shortInt), intent(in) :: binIdx

    call timers(binIdx) % stop()

  end subroutine timerStop

  !!
  !! timerReset Subroutine
  !!
  !! Sets elapsed time to zero
  !!
  !! Args:
  !!   binIdx [in] - Index of the timer, recived from registerTimer function
  !!
  !! Errors:
  !!   If binIdx is out-of-bounds, behaviour is udefined. Most likely SEG error
  !!
  subroutine timerReset(binIdx)
    integer(shortInt), intent(in) :: binIdx

    call timers(binIdx) % reset()

  end subroutine timerReset

  !!
  !! timerTime Sunction
  !!
  !! Returns elapsed time in seconds
  !!
  !! Args:
  !!   binIdx [in] - Index of the timer, recived from registerTimer function
  !!
  !! Result:
  !!   Real with elapsed time in seconds
  !!
  !! Errors:
  !!   If binIdx is out-of-bounds, behaviour is udefined. Most likely SEG error
  !!
  function timerTime(binIdx) result(time)
    integer(shortInt), intent(in) :: binIdx
    real(defReal)                 :: time

    time = timers(binIdx) % toSec()

  end function timerTime

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
