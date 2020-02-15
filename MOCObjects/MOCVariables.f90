module MOCVariables
  use numPrecision
  implicit none

  ! Types for polar quadrature
  ! Presently 3, 4, and 5 aren't supported
  integer(shortInt), parameter, public :: TY = 1, &
                                          GL = 2, &
                                          LO = 3, &
                                          EA = 4, &
                                          EW = 5

end module MOCVariables
