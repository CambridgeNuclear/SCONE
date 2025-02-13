module universalVariables

  use numPrecision

  implicit none

  ! *** DON't CHANGE THIS. HARDCODED IS FINE
  ! CHANGE THIS: NUMBER MUST BE CALCULATED DURING INITIAL GEOMETRY PROCESSING
  ! Problematic for separating modules!
  integer(shortInt), parameter, public :: HARDCODED_MAX_NEST = 8
  integer(shortInt), parameter, public :: MAX_OUTGOING_PARTICLES = 5

  ! CHANGE THIS: NUMBER WILL DEPEND ON SYSTEM ARCHITECTURE
  ! WILL AFFECT PARALLEL SCALING
  integer(shortInt), parameter, public :: array_pad = 64

  ! Display information
  integer(shortInt), parameter, public :: MAX_COL = 70 ! Maximum number of columns in console display

  ! Define variables which are important for tracking neutrons in the geometry
  real(defReal), parameter, public :: INFINITY    = 2.0_defReal**63, &
                                      surface_tol = 1.0e-12_defReal, & ! Tol. on closeness to surface
                                      SURF_TOL    = 1.0E-12_defReal, &
                                      INF         = 2.0_defReal**63, &
                                      NUDGE       = 1.0e-8_defReal     ! Distance to poke neutrons across boundaries for surface tracking

  ! Flags for different possible events in movement in geometry
  integer(shortINt), parameter, public :: COLL_EV = 1, &
                                          BOUNDARY_EV = 2, &
                                          CROSS_EV = 3, &
                                          LOST_EV  = 4

  ! Create definitions for readability when dealing with positions relative to surfaces
  logical(defBool), parameter, public :: behind = .FALSE., &
                                         infront = .TRUE., &
                                         outside = .FALSE., &
                                         inside = .TRUE.

  ! Special material Indexes
  ! NOTE: All material indices MUST BE NON-NEGATIVE!
  integer(shortInt), parameter :: OUTSIDE_MAT = 0 ,&
                                  VOID_MAT    = huge(OUTSIDE_MAT), &
                                  UNDEF_MAT   = VOID_MAT - 1


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
  integer(shortInt), parameter :: X_AXIS = 1 ,&
                                  Y_AXIS = 2 ,&
                                  Z_AXIS = 3

  ! Particle Type Enumeration
  integer(shortInt), parameter :: P_NEUTRON_CE = 1, &
                                  P_NEUTRON_MG = 2

  ! Search error codes
  integer(shortInt), parameter :: valueOutsideArray = -1, &
                                  tooManyIter       = -2, &
                                  targetNotFound    = -3, &
                                  NOT_FOUND         = -3, &
                                  REJECTED          = -4

  ! Integer indexes for type of tracking cross section requested
  integer(shortInt), parameter :: MATERIAL_XS = 1, &
                                  MAJORANT_XS = 2, &
                                  TRACKING_XS = 3

  ! Physical constants
  ! Neutron mass and speed of light in vacuum from from https://physics.nist.gov/cuu/Constants/index.html
  real(defReal), parameter :: neutronMass = 939.56542194_defReal,  & ! Neutron mass in MeV (m*c^2)
                              lightSpeed  = 2.99792458e10_defReal, & ! Light speed in cm/s
                              kBoltzmann  = 1.380649e-23_defReal,  & ! Bolztmann constant in J/K
                              energyPerFission = 200.0_defReal,    &   ! MeV
                              energyDepoZeroScale = 1.04583_defReal ! Ratio of Energy Depos per Fiss in LWR, and Qfiss for U235

  ! Unit conversion
  real(defReal), parameter :: joulesPerMeV = 1.60218e-13  ,&   ! Convert MeV to J
                              shakesPerS   = 1.0e-8            ! Convert shakes to s

  ! Global name variables used to define specific geometry or field types
  character(nameLen), parameter :: nameUFS  = 'uniFissSites'
  character(nameLen), parameter :: nameWW   = 'WeightWindows'

end module universalVariables
