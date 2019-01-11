# Add Source Files to the global list 
add_sources( ./tallyFilter_inter.f90
             ./tallyFilterSlot_class.f90
             ./tallyFilterFactory_func.f90
             ./energyFilter_class.f90
             )
             
add_unit_tests(./Tests/energyFilter_test.f90)