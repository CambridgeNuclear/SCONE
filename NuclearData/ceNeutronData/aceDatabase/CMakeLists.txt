add_sources(./aceNeutronNuclide_class.f90
            ./aceNeutronDatabase_class.f90
            ./urrProbabilityTables_class.f90
	          ./thermalScatteringData_class.f90)

add_integration_tests(./Tests/aceNeutronNuclide_iTest.f90
                      ./Tests/aceNeutronDatabase_iTest.f90
                      ./Tests/urrProbabilityTables_iTest.f90
		                  ./Tests/thermalScatteringData_iTest.f90)

