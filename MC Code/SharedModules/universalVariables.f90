module universalVariables

  use numPrecision

  implicit none

  ! *** DON't CHANGE THIS. HARDCODED IS FINE
  ! CHANGE THIS: NUMBER MUST BE CALCULATED DURING INITIAL GEOMETRY PROCESSING
  ! Problematic for separating modules!
  integer(shortInt), parameter, public :: hardcoded_max_nest = 5

  ! Define variables which are important for tracking neutrons in the geometry
  real(defReal), parameter, public :: INFINITY    = 2.0_defReal**63, &
                                      surface_tol = 1.0e-14_defReal, & ! Tol. on closeness to surface
                                      NUDGE       = 1.0e-8_defReal     ! Distance to poke neutrons across boundaries for surface tracking

  ! Create definitions for readability when dealing with positions relative to surfaces
  logical(defBool), parameter, public :: behind = .FALSE., &
                                         infront = .TRUE., &
                                         outside = .FALSE., &
                                         inside = .TRUE.

  ! Define integers for each fill type that a cell may have
  integer(shortInt), parameter :: OUTSIDE_FILL = 0,  &
                                  materialFill = 1, &
                                  universeFill = 2, &
                                  latticeFill  = 3

  ! Define integers for boundary condition types
  integer(shortInt), parameter :: vacuum     = 0, &
                                  reflective = 1, &
                                  periodic   = 2, &
                                  noBC       = -1
  ! Integer indexes of cardinal directions
  integer(shortInt), parameter :: X_axis = 1 ,&
                                  Y_axis = 2 ,&
                                  Z_axis = 3

  ! Search error codes
  integer(shortInt), parameter :: valueOutsideArray = -1,&
                                  tooManyIter       = -2,&
                                  targetNotFound    = -3

end module universalVariables
