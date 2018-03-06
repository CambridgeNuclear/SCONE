# Define Files to compile 
add_library ( ParticleObjects particle_class.f90
                              particleDungeon_class.f90) 

# Define dependencies of the library 
add_dependencies(ParticleObjects SharedModules Geometry RandomNumbers)