add_subdirectory(keffClerk)
add_subdirectory(tallyMaps)

# Add Source Files to the global list 
add_sources( ./tallyClerk_inter.f90
             ./tallyClerkSlot_class.f90
             ./timeClerk_class.f90
             ./macroClerk_class.f90
             ./tallyClerkFactory_func.f90)