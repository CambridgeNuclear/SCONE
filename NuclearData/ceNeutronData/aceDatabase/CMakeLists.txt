add_sources(./aceNeutronNuclide_class.f90
            ./aceNeutronDatabase_class.f90
            ./aceNeutronNuclideUni_class.f90
            ./aceNeutronDatabaseUni_class.f90
            ./aceNeutronDatabaseUniIdx_class.f90
            ./urrProbabilityTables_class.f90
	          ./thermalScatteringData_class.f90)

add_integration_tests(./Tests/aceNeutronNuclide_iTest.f90
                      ./Tests/aceNeutronDatabase_iTest.f90
                      ./Tests/aceNeutronDatabaseUni_iTest.f90
                      ./Tests/aceNeutronDatabaseUniIdx_iTest.f90
                      ./Tests/urrProbabilityTables_iTest.f90
		                  ./Tests/thermalScatteringData_iTest.f90)

