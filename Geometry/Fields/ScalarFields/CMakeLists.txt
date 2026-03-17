add_subdirectory(PieceConstantFields)

add_sources( ./scalarField_inter.f90
             ./uniformScalarField_class.f90)

add_unit_tests( ./Tests/uniformScalarField_test.f90)
