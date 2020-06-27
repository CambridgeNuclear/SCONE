add_sources(./surface_inter.f90
            ./QuadSurfaces/quadSurface_inter.f90
            ./QuadSurfaces/sphere_class.f90
            ./QuadSurfaces/cylinder_class.f90
            ./QuadSurfaces/aPlane_class.f90
           )

add_unit_tests( ./Tests/sphere_test.f90
                ./Tests/cylinder_test.f90
                ./Tests/aPlane_test.f90
              )
