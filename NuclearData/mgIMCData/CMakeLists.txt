# Add source files for compilation 
add_sources(./mgIMCMaterial_inter.f90
            ./mgIMCDatabase_inter.f90
            ./baseMgIMC/baseMgIMCMaterial_class.f90
            ./baseMgIMC/baseMgIMCDatabase_class.f90) 
            
# Add tests 
#add_integration_tests(./baseMgIMC/Tests/baseMgIMCDatabase_iTest.f90)            
