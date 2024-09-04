module timer_test
  use numPrecision
  use timer_mod, only : registerTimer, timerStart, timerStop, timerReset, timerTime
  use funit

  implicit none


contains

  !!
  !! Test logic for dynamic space for timers
  !! Make sure it does not segment
  !!
@Test
  subroutine testRegisterTimer()
    integer(shortInt) :: i, j

    do i = 1,100
      j = registerTimer('myName')
    end do

    @assertEqual(100, j)

  end subroutine testRegisterTimer

end module timer_test
