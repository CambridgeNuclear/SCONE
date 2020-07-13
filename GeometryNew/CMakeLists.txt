add_subdirectory(Surfaces)
add_subdirectory(Cells)
add_subdirectory(Universes)

add_sources( #./geomSpace_class.f90
             ./geomGraph_class.f90
           )

add_unit_tests( ./Tests/geomGraph_test.f90
              )
