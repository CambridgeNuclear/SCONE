# Add Source Files to the global list 
add_sources(./tallyClerk_inter.f90
            ./tallyClerkFactory_func.f90
            ./tallyClerkSlot_class.f90
            ./collisionClerk_class.f90
            ./keffAnalogClerk_class.f90
            ./keffImplicitClerk_class.f90
            ./simpleFM_class.f90
             )
             
add_unit_tests(./Tests/collisionClerk_test.f90
               ./Tests/keffAnalogClerk_test.f90
               ./Tests/keffImplicitClerk_test.f90
                )