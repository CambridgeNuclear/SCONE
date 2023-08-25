# Add source files for compilation
add_sources(./mgNeutronMaterial_inter.f90
            ./mgNeutronDatabase_inter.f90
            ./baseMgNeutron/baseMgNeutronMaterial_class.f90
            ./baseMgNeutron/baseMgNeutronDatabase_class.f90)

# Add tests
add_integration_tests(./baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90)
