add_subdirectory(muEndfPdf)

# Add Source Files to the global list 
add_sources ( ./angleLawENDF_inter.f90 
              ./angleLawENDFslot_class.f90
              ./angleLawENDFfactory_func.f90
              ./tabularAngle_class.f90
              ./noAngle_class.f90
              ./isotropicAngle_class.f90)



