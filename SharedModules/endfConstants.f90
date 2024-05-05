module endfConstants

  use numPrecision

  implicit none
  integer(shortInt), parameter :: histogramInterpolation = 1 ,&
                                  linLinInterpolation    = 2, &
                                  linLogInterpolation    = 3, &
                                  logLinInterpolation    = 4, &
                                  loglogInterpolation    = 5, &
                                  chargedParticleInterpolation = 6
  ! Seperate parameters for interpolationFlags for tabular PDF. Modern ACE files use flags
  ! consistant with ENDF interpolation paramethers, but old MCNP manual(year 2000) specifies
  ! inconsistant flags: histogram = 0 and linLin = 1
  integer(shortInt), parameter :: tabPdfHistogram = 1 ,&
                                  tabPdfLinLin    = 2

  ! Other constants related to ACE format
  integer(shortInt), parameter :: LOCB_CORRELATED = -1 ,&
                                  LOCB_ISOTROPIC  = 0

  ! List of reaction MT numbers, See Serpent 2 Wiki or ENDF manual for exact details:
  !
  ! Trkov, A., M. Herman, and D. A. Brown. “ENDF-6 Formats Manual.”
  ! Data Formats and Procedures for the Evaluated Nuclear Data Files ENDF/B-VI and ENDF/B-VII,
  ! National Nuclear Data Center Brookhaven National Laboratory, Upton, NY, 2012, 11973–5000.
  !
  ! Dictionary:
  ! N - neutron ; d - deutron; a - alpha particle ; f - fission; p - proton
  !
  ! Syntax:
  ! Example:  N_2N -> (n,2n)
  ! words indicate sum of relevant reactions e.g: (n,fission) is sum (n,f)+(n,nf)+(n,2nf)+(n3nf)

  integer(shortInt), parameter :: N_TOTAL       = 1   ,&
                                  N_N_ELASTIC   = 2   ,&
                                  N_N_INELASTIC = 4   ,&
                                  N_ANYTHING    = 5   ,&
                                  N_2Nd         = 11  ,&
                                  N_2N          = 16  ,&
                                  N_3N          = 17  ,&
                                  N_FISSION     = 18  ,&
                                  N_f           = 19  ,&
                                  N_Nf          = 20  ,&
                                  N_2Nf         = 21  ,&
                                  N_Na          = 22  ,&
                                  N_N3a         = 23  ,&
                                  N_2Na         = 24  ,&
                                  N_3Na         = 25  ,&
                                  N_ABSORPTION  = 27  ,&
                                  N_Np          = 28  ,&
                                  N_N2a         = 29  ,&
                                  N_2N2a        = 30  ,&
                                  N_Nd          = 32  ,&
                                  N_Nt          = 33  ,&
                                  N_Nhe3        = 34  ,&
                                  N_Nd2a        = 35  ,&
                                  N_Nt2a        = 36  ,&
                                  N_4N          = 37  ,&
                                  N_3Nf         = 38  ,&
                                  N_2Np         = 41  ,&
                                  N_3Np         = 42  ,&
                                  N_N2p         = 44  ,&
                                  N_Npa         = 45  ,&
                                  N_Nl1         = 51  ,&
                                 !N_Nl(:) 51-90
                                 !Inelastic scattering from levels 1-40 is defined at the end
                                  N_Nl40        = 90  ,&
                                  N_Ncont       = 91  ,&
                                  N_disap       = 101 ,&
                                  N_GAMMA       = 102 ,&
                                  N_p           = 103 ,&
                                  N_d           = 104 ,&
                                  N_t           = 105 ,&
                                  N_he3         = 106 ,&
                                  N_a           = 107 ,&
                                  N_2a          = 108 ,&
                                  N_3a          = 109 ,&
                                  N_2p          = 111 ,&
                                  N_pa          = 112 ,&
                                  N_t2a         = 113 ,&
                                  N_d2a         = 114 ,&
                                  N_pd          = 115 ,&
                                  N_pt          = 116 ,&
                                  N_da          = 117 ,&
                                  ! SCONE's fake MT for thermal inelastic scattering
                                  N_N_ThermEL     = 1002 ,&
                                  N_N_ThermINEL   = 1004 ,&
                                  ! SCONE's fake MT for particle splitting event
                                  N_N_SPLIT       = 1005

  integer(shortInt),private    :: i  ! Local, private integer to use array constructor
  integer(shortInt),parameter  :: N_Nl(40)      = [(50+i, i =1,40)]
  integer(shortInt),parameter  :: N_2Nl(16)     = [(874+i, i =1,16)]

  ! Microscopic lumped reaction channels special MT numbers
  integer(shortInt),parameter  :: anyScatter    = -102, &
                                  anyCapture    = -201, &
                                  anyFission    = -118

  ! List of Fake MT numbers for macroscopic XSs. Stolen from Serpent
  integer(shortInt),parameter  :: macroTotal     = -1 ,&
                                  macroCapture   = -2 ,&
                                  macroEscatter  = -3 ,&
                                  macroIEscatter = -4 ,&
                                  macroFission   = -6 ,&
                                  macroNuFission = -7

  ! List of Macro MT numbers for macroscopic XSs. Unique to SCONE (not from Serpent)
  integer(shortInt), parameter :: macroAllScatter = -20 ,&
                                  macroAbsorbtion = -21 ,&
                                  noInteraction   = -901






  ! List of diffrent ENDF energy Laws. See ENDF manual for more details:
  !
  ! Trkov, A., M. Herman, and D. A. Brown. “ENDF-6 Formats Manual.”
  ! Data Formats and Procedures for the Evaluated Nuclear Data Files ENDF/B-VI and ENDF/B-VII,
  ! National Nuclear Data Center Brookhaven National Laboratory, Upton, NY, 2012, 11973–5000.
  !
  integer(shortInt),parameter  :: tabularEquiprobableEnergyBins = 1  ,&
                                  discretePhotonEnergy          = 2  ,&
                                  levelScatteringLaw            = 3  ,&
                                  continuousTabularDistribution = 4  ,&
                                  generalEvaporationSpectrum    = 5  ,&
                                  simpleMaxwellFissionSpectrum  = 7  ,&
                                  evaporationEnergySpectrum     = 9  ,&
                                  energyDependentWattSpectrum   = 11 ,&
                                  tabularLinearFunctions        = 22 ,&
                                  ukLaw6                        = 24 ,&
                                  kalbach87Formalism            = 44 ,&
                                  endfEnergyLaw61               = 61 ,&
                                  nBodyPhaseSpaceDistribution   = 66 ,&
                                  labratoryAngleEnergyLaw       = 67


end module endfConstants
