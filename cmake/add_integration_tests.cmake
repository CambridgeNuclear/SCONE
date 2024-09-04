# Accumulate source files in a global property INTEGRATION_TEST_LIST
# Store source files with their absolute files
#
# At the end INTEGRATION_TEST_LIST can be converted to a variable and used to
# compile large number of files spread across multiple folders with a single command
#
function(add_integration_tests)
  # Check if the property is already defined
  # If it isn't define it
  get_property(is_defined GLOBAL PROPERTY INTEGRATION_TEST_LIST DEFINED)
  if(NOT is_defined)
    define_property(GLOBAL PROPERTY INTEGRATION_TEST_LIST
      BRIEF_DOCS "List of all integration test files"
      FULL_DOCS "List of all pFUnit integration test suite files for preprocessing & compilation")
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
  set_property(GLOBAL APPEND PROPERTY INTEGRATION_TESTS_LIST "${TESTS}")

endfunction(add_integration_tests)
