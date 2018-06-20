add_subdirectory(Planes)
add_subdirectory(Cylinders)
add_subdirectory(SquareCylinders)
add_subdirectory(TruncatedCylinders) 

# Add Source Files to the global list 
add_sources( ./box_class.f90
             ./infSurf_class.f90
             ./sphere_class.f90
             ./surface_inter.f90
             ./surfaceFactory_func.f90) 