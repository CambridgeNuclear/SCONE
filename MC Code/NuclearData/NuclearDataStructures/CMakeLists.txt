# Define Files to compile 
add_library( NuclearDataStructures ./endfTable/endfTable_class.f90
                                   ./pdf/tabularPdf_class.f90)

# Define dependencies of the library 
target_link_libraries(NuclearDataStructures SharedModules )