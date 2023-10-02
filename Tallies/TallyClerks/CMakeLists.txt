# Add Source Files to the global list 
add_sources(./tallyClerk_inter.f90
            ./tallyClerkFactory_func.f90
            ./tallyClerkSlot_class.f90
            ./collisionClerk_class.f90
            ./collisionProbabilityClerk_class.f90
            ./trackClerk_class.f90
            ./keffAnalogClerk_class.f90
            ./keffImplicitClerk_class.f90
            ./simpleFMClerk_class.f90
            ./dancoffBellClerk_class.f90
            ./shannonEntropyClerk_class.f90
            ./centreOfMassClerk_class.f90
            ./mgXsClerk_class.f90
             )
             
add_unit_tests(./Tests/collisionClerk_test.f90
               ./Tests/trackClerk_test.f90
               ./Tests/keffAnalogClerk_test.f90
               ./Tests/keffImplicitClerk_test.f90
               ./Tests/shannonEntropyClerk_test.f90
               ./Tests/simpleFMClerk_test.f90
               ./Tests/collisionProbabilityClerk_test.f90
               ./Tests/mgXsClerk_test.f90
                )
