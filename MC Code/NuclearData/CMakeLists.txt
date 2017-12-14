# Define Files to compile 
add_library( NuclearData nuclearData_class.f90 
                         ./ByIsoNoMT/ByIsoNoMT_Data_class.f90 
                         ./ByIsoNoMT/aceNoMT_class.f90)

# Define dependencies of the library 
target_link_libraries(NuclearData SharedModules RandomNumbers ParticleObjects)