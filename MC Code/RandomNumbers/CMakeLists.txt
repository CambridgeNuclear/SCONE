# Define Files to compile 
add_library(RandomNumbers RNG_class.f90)

# Define dependencies of the library 
target_link_libraries(RandomNumbers SharedModules)