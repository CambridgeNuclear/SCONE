# Add Nested Directories 
add_subdirectory(angleLawENDF)
#add_subdirectory(energyLawENDF)

# Add Source Files to the global list 
add_sources( ./emissionENDF_class.f90
             ./correlatedEmissionENDF_class.f90 
             ./uncorrelatedEmissionENDF_class.f90
             ./emissionFromACE_func.f90
             ./releaseLawENDF/releaseLawENDF_class.f90
             ./releaseLawENDF/polynomialRelease_class.f90
             ./releaseLawENDF/constantRelease_class.f90
             ./releaseLawENDF/tabularRelease_class.f90
             ./energyLawENDF/energyLawENDF_class.f90
             ./energyLawENDF/contTabularEnergy_class.f90
             ./energyLawENDF/energyEndfPdf/tabularEnergy_class.f90
             ./energyLawENDF/noEnergy_class.f90
             ./energyLawENDF/maxwellSpectrum_class.f90)
