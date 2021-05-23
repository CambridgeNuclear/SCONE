# Add Source Files to the global list
add_sources( ./outputFile_class.f90
             ./asciiOutput_inter.f90
             ./asciiOutputFactory_func.f90
             ./asciiMATLAB_class.f90
             ./asciiJSON_class.f90
             ./dummyPrinter_class.f90)

add_unit_tests (./Tests/outputFile_test.f90)
