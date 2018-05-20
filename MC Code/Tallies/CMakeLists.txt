#add_subdirectory(TallyResponses)
#add_subdirectory(TallyProc)

add_subdirectory(TallyClerks)

# Add Source Files to the global list 
add_sources( ./tallyAdminBase_class.f90
             ./tallyCounterOLD_class.f90
             ./tallyEstimator_class.f90 )
             
