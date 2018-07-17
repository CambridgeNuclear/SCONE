add_subdirectory(Planes)
add_subdirectory(Cylinders)
#add_subdirectory(TruncCylinders) 

# Add Source Files to the global list 
add_sources( ./surface_inter.f90
             ./sphere_class.f90 
             #./surfaceFactory_func.f90
             ./box_class.f90
             ./squareCylinder_class.f90
             ./infSurf_class.f90
             #./sphere_class.f90) 
             )