# Add Nested Directories 
add_subdirectory(angleLawENDF)
add_subdirectory(energyLawENDF)
add_subdirectory(releaseLawENDF) 
add_subdirectory(correlatedLawENDF) 

# Add Source Files to the global list 
add_sources( ./emissionENDF_class.f90) 
