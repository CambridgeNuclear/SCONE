#add_subdirectory(TallyResponses)
#add_subdirectory(TallyProc)

add_subdirectory(TallyClercs)

# Add Source Files to the global list 
add_sources( ./tallyAdminBase_class.f90
             ./tallyCounter_class.f90 )
             
