add_sources( ./cell_inter.f90
             ./cellFactory_func.f90
             ./cellShelf_class.f90
             ./simpleCell_class.f90
           )

add_unit_tests( ./Tests/simpleCell_test.f90
                ./Tests/cellShelf_test.f90
              )
