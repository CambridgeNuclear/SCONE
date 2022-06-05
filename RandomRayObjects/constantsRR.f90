module constantsRR
  
  use numPrecision
  
  implicit none

  ! Parameter for when to skip a tiny volume
  real(defReal), parameter, public :: volume_tolerance = 1.0E-12
  
  ! Parameters to identify the simulation type
  integer(shortInt), parameter, public :: flatIso = 1, linearIso = 2, flatAni = 3, linearAni = 4


contains
    
end module constantsRR
