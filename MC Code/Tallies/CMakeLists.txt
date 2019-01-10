#add_subdirectory(TallyResponses)
#add_subdirectory(TallyProc)

add_subdirectory(TallyClerks)
add_subdirectory(TallyFilters)
add_subdirectory(TallyMaps)

# Add Source Files to the global list 
add_sources( ./tallyCodes.f90
             ./tallyAdminBase_class.f90
             #./tallyInactiveAdmin_class.f90
             #./tallyActiveAdmin_class.f90 
             #./tallyTimeAdmin_class.f90
             ./tallyEstimator_class.f90 
             ./scoreMemory_class.f90)
             
add_unit_tests(./Tests/scoreMemory_test.f90)
                     
             
