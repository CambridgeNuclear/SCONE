# Add Nested Directories 
add_subdirectory(angleLawENDF)
add_subdirectory(energyLawENDF)
add_subdirectory(releaseLawENDF) 

# Add Source Files to the global list 
add_sources( ./emissionENDF_class.f90
             ./correlatedEmissionENDF_class.f90 
             ./uncorrelatedEmissionENDF_class.f90
             ./emissionFromACE_func.f90)
