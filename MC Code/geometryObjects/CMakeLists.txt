# Define Files to compile 
add_library ( Geometry coord_class.f90) 

# Define dependencies of the library 
add_dependencies(Geometry SharedModules )