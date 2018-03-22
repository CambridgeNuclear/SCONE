# Define Files to compile 
add_library ( CollisionOperator ./collisionOperator_class.f90
                                ./scatteringKernels_func.f90
                                ./MGOperators/collisionOperatorMG_class.f90) 

# Define dependencies of the library 
add_dependencies(CollisionOperator SharedModules NuclearData RandomNumbers ParticleObjects )