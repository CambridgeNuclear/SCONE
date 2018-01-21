# Define Files to compile 
add_library( angleLawENDF angleLawENDF_class.f90
                          tabularAngle_class.f90
                          noAngle_class.f90
                          ./miuEndfPdf/miuEndfPdf_class.f90
                          ./miuEndfPdf/isotropicMiu_class.f90
                          ./miuEndfPdf/equiBin32Miu_class.f90
                          ./miuEndfPdf/tabularMiu_class.f90)

 

# Define dependencies of the library 
target_link_libraries(angleLawENDF SharedModules RandomNumbers NuclearDataStructures )