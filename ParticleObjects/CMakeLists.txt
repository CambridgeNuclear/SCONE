add_subdirectory(Source)

# Add Source Files to the global list 
add_sources( particle_class.f90
             particleDungeon_class.f90 ) 

add_unit_tests( ./Tests/particle_test.f90
                ./Tests/particleDungeon_test.f90) 
