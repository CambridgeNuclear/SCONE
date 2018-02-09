add_subdirectory(angleLawENDF)
#add_subdirectory(energyLawENDF)
# Define Files to compile 
add_library( emissionENDF emissionENDF_class.f90
                          correlatedEmissionENDF_class.f90 
                          uncorrelatedEmissionENDF_class.f90
                          emissionFromACE_func.f90
                          ./releaseLawENDF/releaseLawENDF_class.f90
                          ./releaseLawENDF/polynomialRelease_class.f90
                          ./releaseLawENDF/constantRelease_class.f90
                          ./releaseLawENDF/tabularRelease_class.f90
                          ./energyLawENDF/energyLawENDF_class.f90
                          ./energyLawENDF/contTabularEnergy_class.f90
                          ./energyLawENDF/energyEndfPdf/tabularEnergy_class.f90
                          ./energyLawENDF/noEnergy_class.f90
                          ./energyLawENDF/maxwellSpectrum_class.f90)

 

# Define dependencies of the library 
target_link_libraries(emissionENDF RandomNumbers ParticleObjects NuclearDataStructures 
                                   angleLawENDF )