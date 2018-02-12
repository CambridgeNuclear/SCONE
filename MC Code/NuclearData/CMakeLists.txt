add_subdirectory(emissionENDF)
add_subdirectory(NuclearDataStructures) 
# Define Files to compile 
add_library( NuclearData nuclearData_class.f90 
                         ./ByIsoNoMT/ByIsoNoMT_Data_class.f90 
                         ./ByIsoNoMT/aceNoMT_class.f90
                         ./ByIsoNoMT/ByIsoNoMT_class.f90
                         ./CrossSectionPackages/xsCDF_class.f90
                         ./CrossSectionPackages/xsMainCDF_class.f90 
                         ./CrossSectionPackages/xsAnyCDF_class.f90
                         ./CrossSectionPackages/xsSet_class.f90)
# Define dependencies of the library 
target_link_libraries(NuclearData SharedModules RandomNumbers ParticleObjects emissionENDF)