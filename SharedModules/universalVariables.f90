module universalVariables

  use numPrecision

  implicit none

  ! *** DON't CHANGE THIS. HARDCODED IS FINE
  ! CHANGE THIS: NUMBER MUST BE CALCULATED DURING INITIAL GEOMETRY PROCESSING
  ! Problematic for separating modules!
  integer(shortInt), parameter, public :: hardcoded_max_nest = 5
  integer(shortInt), parameter, public :: MAX_OUTGOING_PARTICLES = 5

  ! Define variables which are important for tracking neutrons in the geometry
  real(defReal), parameter, public :: INFINITY    = 2.0_defReal**63, &
                                      surface_tol = 1.0e-12_defReal, & ! Tol. on closeness to surface
                                      SURF_TOL    = 1.0E-12_defReal, &
                                      INF         = 2.0_defReal**63, &
                                      NUDGE       = 1.0e-8_defReal     ! Distance to poke neutrons across boundaries for surface tracking

  ! Create definitions for readability when dealing with positions relative to surfaces
  logical(defBool), parameter, public :: behind = .FALSE., &
                                         infront = .TRUE., &
                                         outside = .FALSE., &
                                         inside = .TRUE.

  ! Special material Indexes
  ! NOTE: All material indices MUST BE NON-NEGATIVE!
  integer(shortInt), parameter :: OUTSIDE_MAT = 0 ,&
                                  VOID_MAT    = huge(OUTSIDE_MAT)


  ! Define integers for each fill type that a cell may have
  integer(shortInt), parameter :: OUTSIDE_FILL = 0,  &
                                  materialFill = 1, &
                                  universeFill = 2, &
                                  latticeFill  = 3

  ! Define integers for boundary condition types
  integer(shortInt), parameter :: VACUUM_BC     = 0, &
                                  REFLECTIVE_BC = 1, &
                                  PERIODIC_BC   = 2
                                  
  ! Integer indexes of cardinal directions
  integer(shortInt), parameter :: X_axis = 1 ,&
                                  Y_axis = 2 ,&
                                  Z_axis = 3

  ! Particle Type Enumeration
  integer(shortInt), parameter :: P_NEUTRON_CE = 1, &
                                  P_NEUTRON_MG = 2

  ! Search error codes
  integer(shortInt), parameter :: valueOutsideArray = -1,&
                                  tooManyIter       = -2,&
                                  targetNotFound    = -3, &
                                  NOT_FOUND         = -3

  ! Physical constants
  real(defReal), parameter :: neutronMass = 939.5654133_defReal, &   ! Neutron mass in MeV/c^2
                              lightSpeed  = 2.99792458e10_defReal, & ! Light speed in cm/s
                              energyPerFission = 200.0_defReal       ! MeV

  ! Unit conversion
  real(defReal), parameter :: joulesPerMeV = 1.60218e-13     ! Convert MeV to J

end module universalVariables
