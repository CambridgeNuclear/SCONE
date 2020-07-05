add_sources( ./universe_inter.f90
             ./cellUniverse_class.f90
             ./pinUniverse_class.f90
           )


add_unit_tests( ./Tests/universe_test.f90
                ./Tests/cellUniverse_test.f90
                ./Tests/pinUniverse_test.f90
              )
