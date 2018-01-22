# Define Files to compile 
add_library( angleLawENDF angleLawENDF_class.f90
                          tabularAngle_class.f90
                          noAngle_class.f90
                          ./muEndfPdf/muEndfPdf_class.f90
                          ./muEndfPdf/isotropicMu_class.f90
                          ./muEndfPdf/equiBin32Mu_class.f90
                          ./muEndfPdf/tabularMu_class.f90)

 

# Define dependencies of the library 
target_link_libraries(angleLawENDF SharedModules RandomNumbers NuclearDataStructures )