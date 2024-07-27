module constantsRR
  
  use numPrecision
  
  implicit none
  
  ! Parameters for volume/no-hit policy
  integer(shortInt), parameter :: hybrid = 3, srcPolicy = 1, prevPolicy = 2, &
                                  naive = 2, simAverage = 1

  ! Parameters to identify the simulation type
  integer(shortInt), parameter, public :: flatIso = 1, linearIso = 2, flatAni = 3, linearAni = 4
  
  ! Parameter for when to skip a tiny volume
  real(defReal), parameter, public :: volume_tolerance = 1.0E-12

  ! Parameter for when to ignore components of spatial moment matrices
  ! or when the matrix is poorly conditioned
  real(defReal), parameter, public :: condition_tolerance = 1.0E-10, &
                                      det_tolerance = 1.0E-12
  
  ! Parameters for indexing into matrices and spatial moments with linear sources
  integer(shortInt), parameter :: x = 1, y = 2, z = 3, nDim = 3, &
                                  xx = 1, xy = 2, xz = 3, &
                                  yy = 4, yz = 5, zz = 6, &
                                  matSize = 6

  ! Parameters for deciding how to invert the moment matrix
  integer(shortInt), parameter :: invertXYZ = 7, invertXY = 6, &
                                  invertXZ = 5, invertYZ = 3, &
                                  invertX = 4, invertY = 2, &
                                  invertZ = 1

  ! Convenient arithmetic parameters
  real(defFlt), parameter :: one_two = real(HALF,defFlt), &
                             two_three = real(2.0_defFlt/3.0_defFlt,defFlt)

  real(defReal), parameter :: one_twelve = ONE / 12


contains
    
end module constantsRR
