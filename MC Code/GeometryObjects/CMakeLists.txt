add_subdirectory(SurfaceObjects)

# Add Source Files to the global list 
add_sources( ./geometry_inter.f90 
             ./cellGeometry_inter.f90
             ./genGeometry_inter.f90
            # ./cell_class.f90
             ./coord_class.f90)
            # ./geometry_class.f90
            # ./lattice_class.f90
            # ./mesh_class.f90
            # ./universe_class.f90)
