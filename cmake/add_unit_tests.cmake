# Accumulate source files in a global property UNIT_TEST_LIST
# Store source files with their absolute files
#
# At the end UNIT_TEST_LIST can be converted to a variable and used to
# compile large number of files spread across multiple folders with a single command
#
function(add_unit_tests)
  # Check if the property is already defined
  # If it isn't define it
  get_property(is_defined GLOBAL PROPERTY UNIT_TEST_LIST DEFINED)
  if(NOT is_defined)
    define_property(GLOBAL PROPERTY UNIT_TEST_LIST
      BRIEF_DOCS "List of all test files"
      FULL_DOCS "List of all pFUnit test suite files for preprocessing & compilation")
  endif()

  # Take files listed in argument list and make their paths absolute
  set(TESTS)
  foreach(t IN LISTS ARGN)
    if(NOT IS_ABSOLUTE "${t}")
      get_filename_component(t "${t}" ABSOLUTE)
    endif()
    list(APPEND TESTS "${t}")
  endforeach()

  # Append files in argument list to global property
  set_property(GLOBAL APPEND PROPERTY UNIT_TESTS_LIST "${TESTS}")

endfunction(add_unit_tests)
