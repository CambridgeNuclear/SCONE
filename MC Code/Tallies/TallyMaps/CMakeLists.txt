# Add Source Files to the global list 
add_sources(./tallyMap_inter.f90
            ./tallyMap1D_inter.f90 
 #           ./tallyMapFactory_func.f90
 #           ./tallyMapSlot_class.f90
 #           ./energyMap_class.f90
 #           ./spaceMap_class.f90
            ./materialMap_class.f90
 #           ./matXsMap_class.f90
            )
            
add_unit_tests(./Tests/materialMap_test.f90) 