add_subdirectory(ceNeutronData)
add_subdirectory(mgNeutronData)
add_subdirectory(testNeutronData)
add_subdirectory(xsPackages)
add_subdirectory(emissionENDF)
add_subdirectory(NuclearDataStructures) 

add_sources(./nuclearDatabase_inter.f90 
            ./materialHandle_inter.f90
            ./nuclideHandle_inter.f90
            ./materialMenu_mod.f90
            ./nuclearDataReg_mod.f90
            ./neutronMaterial_inter.f90
            ./Reactions/reactionHandle_inter.f90
            ./Reactions/uncorrelatedReactionCE_inter.f90
            ./Reactions/reactionMG_inter.f90
            ./Reactions/correlatedReactionCE_inter.f90
            ./Reactions/uncorrelatedReactionCE/elasticNeutronScatter_class.f90
            ./Reactions/uncorrelatedReactionCE/neutronScatter_class.f90 
            ./Reactions/uncorrelatedReactionCE/fissionCE_class.f90 
            ./Reactions/uncorrelatedReactionCE/pureCapture_class.f90
            ./Reactions/reactionMG/fissionMG_class.f90
            ./Reactions/reactionMG/multiScatterMG_class.f90 
            ./Reactions/reactionMG/multiScatterP1MG_class.f90 
            ./DataDecks/dataDeck_inter.f90
            ./DataDecks/ACE/aceCard_class.f90
            ./DataDecks/dictDeck_class.f90) 

add_unit_tests(./Reactions/reactionMG/Tests/fissionMG_test.f90
               ./Reactions/reactionMG/Tests/multiScatterMG_test.f90
               ./Reactions/reactionMG/Tests/multiScatterP1MG_test.f90
               ./Reactions/Tests/pureCapture_test.f90
               ./Tests/materialMenu_test.f90) 

            
add_integration_tests(./Reactions/Tests/elasticScattering_iTest.f90
                      ./Reactions/Tests/neutronScattering_iTest.f90
                      ./Reactions/Tests/fissionCE_iTest.f90
                      ./Tests/aceLibraryRead_iTest.f90) 
