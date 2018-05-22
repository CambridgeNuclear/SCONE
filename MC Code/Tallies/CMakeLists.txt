#add_subdirectory(TallyResponses)
#add_subdirectory(TallyProc)

add_subdirectory(TallyClerks)

# Add Source Files to the global list 
add_sources( ./tallyCodes.f90
             ./tallyAdminBase_class.f90
             ./tallyInactiveAdmin_class.f90 
             ./tallyEstimator_class.f90 )
             
