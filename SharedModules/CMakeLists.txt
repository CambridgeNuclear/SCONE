# Add Source Files to the global list
add_sources( ./genericProcedures.f90
             ./numPrecision.f90
             ./endfConstants.f90
             ./universalVariables.f90
             ./hashFunctions_func.f90
             ./legendrePoly_func.f90
             ./charTape_class.f90
             ./grid_class.f90
             ./energyGrid_class.f90
             ./statisticalTests_func.f90
	     ./exponentialRA_func.f90
             ./timer_mod.f90
             ./charLib_func.f90
             ./openmp_func.f90
             ./errors_mod.f90
             ./colours_func.f90)

add_unit_tests( ./Tests/grid_test.f90
                ./Tests/energyGrid_test.f90
                ./Tests/sort_test.f90
                ./Tests/statisticalTests_test.f90
		./Tests/exponentialRA_test.f90
                ./Tests/hashFunctions_test.f90
                ./Tests/timer_test.f90
                ./Tests/conversions_test.f90
                ./Tests/charLib_test.f90
                ./Tests/colours_test.f90)
