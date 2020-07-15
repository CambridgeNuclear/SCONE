add_subdirectory(Surfaces)
add_subdirectory(Cells)
add_subdirectory(Universes)
add_subdirectory(Fields)

add_sources( ./csg_class.f90
             ./geomGraph_class.f90
             ./geometry_inter.f90
             ./geometryStd_class.f90
           )

add_unit_tests( ./Tests/geomGraph_test.f90
              )

add_integration_tests( ./Tests/geometryStd_iTest.f90
                     )
