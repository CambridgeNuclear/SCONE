module numPrecision
  implicit none
  private
  ! Variables Kind and Length parameters
  integer, public, parameter :: defReal = 8,     &
                                defFlt  = 8,     &
                                shortInt = 4,    &
                                longInt = 8,     &
                                defBool = 4,     &
                                pathLen = 100,   &
                                nameLen = 30
  ! I/O error codes
  integer, public, parameter :: endOfFile = -1


  ! Usefull constants
  real(defReal), public, parameter :: PI = 4.0_defReal * atan(1.0_defReal), &
                                      SQRT2 = sqrt(2._defReal), &
                                      SQRT2_2 = sqrt(2._defReal)/2._defReal , &
                                      ZERO = 0._defReal, &
                                      ONE = 1.0_defReal, &
                                      TWO = 2.0_defReal, &
                                      TWO_PI = TWO * PI, &
                                      HALF = 0.5_defReal,&
                                      FOUR_PI = TWO * TWO_PI, &
                                      ONE_FOUR_PI = ONE /(FOUR_PI)

  real(defReal), public, parameter  :: floatTol = 1.0e-12 !*** Should be replaced
  real(defReal), public, parameter  :: FP_REL_TOL = 1.0e-7_defReal


contains
    
end module numPrecision
