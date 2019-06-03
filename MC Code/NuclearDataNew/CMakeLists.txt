add_sources(./Reactions/reactionHandle_inter.f90
            ./Reactions/uncorrelatedReactionCE_inter.f90
            ./Reactions/reactionMG_inter.f90
            ./Reactions/correlatedReactionCE_inter.f90
            ./Reactions/uncorrelatedReactionCE/elasticNeutronScatter_class.f90
            ./Reactions/uncorrelatedReactionCE/neutronScatter_class.f90 
            ./DataDecks/dataDeck_inter.f90
            ./DataDecks/ACE/aceCard_class.f90) 
            
add_integration_tests(./Reactions/Tests/elasticScattering_iTest.f90
                      ./Reactions/Tests/neutronScattering_iTest.f90) 