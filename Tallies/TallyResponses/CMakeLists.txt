# Add Source Files to the global list 
add_sources(./tallyResponse_inter.f90
            ./tallyResponseFactory_func.f90 
            ./tallyResponseSlot_class.f90
            ./fluxResponse_class.f90
            ./macroResponse_class.f90
            ./weightResponse_class.f90
            ./testResponse_class.f90
            )

             
add_unit_tests(./Tests/fluxResponse_test.f90
               ./Tests/testResponse_test.f90 
               ./Tests/macroResponse_test.f90
               ./Tests/weightResponse_test.f90
                )