module numPrecision
  implicit none
  private
  ! Variables Kind and Length parameters
  integer, public, parameter :: defReal = 8,     &
                                shortInt = 4,    &
                                longInt = 8,     &
                                defBool = 4,     &
                                matNameLen = 20, &
                                pathLen = 100,   &
                                zzIdLen = 10,    &
                                lineLen = 300,   &
                                nameLen = 30
  ! I/O error codes
  integer, public, parameter :: endOfFile = -1


  ! Usefull constants
  real(defReal), public, parameter :: PI = 4.0 * atan(1.0_defReal), &
                                      SQRT2 =sqrt(2._defReal), &
                                      SQRT2_2 = sqrt(2._defReal)/2._defReal , &
                                      ZERO = 0._defReal, &
                                      ONE = 1.0_defReal, & ! CHeck wheter ONE or TWO are used
                                      TWO = 2.0_defReal

  !real(defReal), public, parmater  :: floatTol = 1.0e-23

contains
    
end module numPrecision
