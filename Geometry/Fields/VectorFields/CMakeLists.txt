
add_sources( ./vectorField_inter.f90
             ./uniformVectorField_class.f90
             ./weightWindowsField_class.f90)

add_unit_tests( ./Tests/uniformVectorField_test.f90
                ./Tests/weightWindowsField_test.f90)
