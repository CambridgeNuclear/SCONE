add_subdirectory(emissionENDF)
add_subdirectory(NuclearDataStructures) 
# Define Files to compile 
add_library( NuclearData nuclearData_class.f90 
                         ./byNucNoMT/byNucNoMT_Data_class.f90 
                         ./byNucNoMT/aceNoMT_class.f90
                         ./byNucNoMT/byNucNoMT_class.f90
                         ./byNucNoMT/xsEnergyPointNoMT_class.f90
                         ./byNucNoMT/nuclideMemoryNoMT_class.f90   
                         ./byNucNoMT/materialDataNoMT_class.f90  
                         ./byNucNoMT/materialMemoryNoMT_class.f90                         
                         ./CrossSectionPackages/xsCDF_class.f90
                         ./CrossSectionPackages/xsMainCDF_class.f90 
                         ./CrossSectionPackages/xsAnyCDF_class.f90
                         ./CrossSectionPackages/xsSet_class.f90
                         ./CrossSectionPackages/xsMainSet_class.f90
                         ./CrossSectionPackages/matNucCDF_class.f90                         
                         ./CrossSectionPackages/xsMacroSet_class.f90
                         ./MGnuclearData/perMaterialMgXs_inter.f90)
# Define dependencies of the library 
target_link_libraries(NuclearData SharedModules RandomNumbers ParticleObjects emissionENDF DataStructures)