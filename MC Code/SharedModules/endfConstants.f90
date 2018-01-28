module endfConstants
  implicit none
  integer, parameter :: histogramInterpolation = 1 ,&
                        linLinInterpolation = 2, &
                        linLogInterpolation = 3, &
                        logLinInterpolation = 4, &
                        loglogInterpolation = 5, &
                        chargedParticleInterpolation = 6
  ! Seperate parameters for interpolationFlags for tabular PDF. Modern ACE files use flags
  ! consistant with ENDF interpolation paramethers, but old MCNP manual(year 2000) specifies
  ! inconsistant flags: histogram = 0 and linLin = 1
  integer, parameter :: tabPdfHistogram = 1 ,&
                        tabPdfLinLin    = 2

end module endfConstants
