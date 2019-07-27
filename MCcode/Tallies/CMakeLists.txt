add_subdirectory(TallyClerks)
add_subdirectory(TallyFilters)
add_subdirectory(TallyMaps)
add_subdirectory(TallyResponses) 

# Add Source Files to the global list 
add_sources( ./tallyCodes.f90
             ./tallyResult_class.f90
             ./tallyAdmin_class.f90
             ./scoreMemory_class.f90)
             
add_unit_tests(./Tests/scoreMemory_test.f90)
                     
             
