# Add Source Files to a global List
add_sources( ./dictionary_class.f90
             ./intMap_class.f90
             ./charMap_class.f90
             ./stack_class.f90
             ./dynArray_class.f90
             ./dictParser_func.f90
             ./heapQueue_class.f90
             ./linkedList_class.f90)

# Add Unit Tests to a global List
add_unit_tests( ./Tests/dictionary_test.f90
                ./Tests/intMap_test.f90
                ./Tests/charMap_test.f90
                ./Tests/dynArray_test.f90
                ./Tests/dictParser_test.f90
                ./Tests/heapQueue_test.f90
                ./Tests/linkedList_test.f90)

add_integration_tests( ./Tests/dictParser_iTest.f90)
