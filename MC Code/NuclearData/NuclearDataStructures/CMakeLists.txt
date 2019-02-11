add_sources ( ./endfTable/endfTable_class.f90
              ./pdf/tabularPdf_class.f90
              ./pdf/maxwellEnergyPdf_class.f90
              ./pdf/kalbachTable_class.f90)

add_unit_tests(./Tests/tabularPdf_test.f90)