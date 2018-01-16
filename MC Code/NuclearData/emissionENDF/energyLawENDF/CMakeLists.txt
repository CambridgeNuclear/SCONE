# Define Files to compile 
add_library( energyLawEndf energyLawENDF_class.f90)
                         

 
# Define dependencies of the library 
target_link_libraries(energyLawEndf SharedModules RandomNumbers NuclearDataStructures)