# Add Nested Directories 
add_subdirectory(emissionENDF)
add_subdirectory(NuclearDataStructures) 
add_subdirectory(transportNuclearData) 
add_subdirectory(AceParser) 

# Add Source Files to the global list 
add_sources( ./nuclearData_inter.f90
             ./nuclearDataRegistry_mod.f90
             ./datalessMaterials_class.f90                     
             ./CrossSectionPackages/xsCDF_class.f90
             ./CrossSectionPackages/xsMainCDF_class.f90 
             ./CrossSectionPackages/xsAnyCDF_class.f90
             ./CrossSectionPackages/xsMainSet_class.f90
             ./CrossSectionPackages/xsNucMacroSet_class.f90                         
             ./CrossSectionPackages/xsMacroSet_class.f90)
