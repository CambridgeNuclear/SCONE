add_sources(./surface_inter.f90
            ./QuadSurfaces/quadSurface_inter.f90
            ./QuadSurfaces/sphere_class.f90
            ./QuadSurfaces/cylinder_class.f90
            ./QuadSurfaces/aPlane_class.f90
            ./QuadSurfaces/plane_class.f90
            ./QuadSurfaces/cone_class.f90
            ./QuadSurfaces/quadric_class.f90
            ./box_class.f90
            ./squareCylinder_class.f90
            ./truncCylinder_class.f90
            ./surfaceFactory_func.f90
            ./surfaceShelf_class.f90
           )

add_unit_tests( ./Tests/sphere_test.f90
                ./Tests/cylinder_test.f90
                ./Tests/cone_test.f90
                ./Tests/aPlane_test.f90
                ./Tests/plane_test.f90
                ./Tests/quadric_test.f90
                ./Tests/box_test.f90
                ./Tests/squareCylinder_test.f90
                ./Tests/truncCylinder_test.f90
                ./Tests/surfaceShelf_test.f90
              )
