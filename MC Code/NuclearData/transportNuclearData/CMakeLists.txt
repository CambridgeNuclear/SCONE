# Add Nested Directories 
add_subdirectory(byNucNoMT)
add_subdirectory(byNucMT)
add_subdirectory(perMaterialNuclearDataMG)


# Add Source Files to the global list 
add_sources( ./transportNuclearData_inter.f90
             ./perNuclideNuclearDataCE_inter.f90
             ./perMaterialNuclearDataMG_inter.f90) 