# Add Source Files to the global list
add_sources(./tallyMap_inter.f90
            ./tallyMap1D_inter.f90
            ./tallyMap1DFactory_func.f90
            ./tallyMapFactory_func.f90
            ./tallyMapSlot_class.f90
            ./testMap_class.f90
            ./energyMap_class.f90
            ./spaceMap_class.f90
            ./materialMap_class.f90
            ./homogMatMap_class.f90
            ./multiMap_class.f90
            ./weightMap_class.f90
            ./sphericalMap_class.f90
            ./cellMap_class.f90
            ./cylindricalMap_class.f90
 #          ./matXsMap_class.f90
            )

add_unit_tests(./Tests/materialMap_test.f90
               ./Tests/energyMap_test.f90
               ./Tests/weightMap_test.f90
               ./Tests/spaceMap_test.f90
               ./Tests/testMap_test.f90
               ./Tests/multiMap_test.f90
               ./Tests/homogMatMap_test.f90
               ./Tests/sphericalMap_test.f90
               ./Tests/cellMap_test.f90)
