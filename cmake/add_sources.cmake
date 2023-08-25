# Accumulate source files in a global property SRCS_LIST
# Store source files with their absolute files
#
# Add the end SCRS_LIST can be converted to a variable and used to
# compile large number of files spread across multiple folders with a single command
#
function(add_sources)
  # Check if the property is already defined
  # If it isn't define it
  get_property(is_defined GLOBAL PROPERTY SRCS_LIST DEFINED)
  if(NOT is_defined)
    define_property(GLOBAL PROPERTY SRCS_LIST
      BRIEF_DOCS "List of source files"
      FULL_DOCS "List of source files to be compiled in one library")
  endif()

  # Take files listed in argument list and make their paths absolute
  set(SRCS)
  foreach(s IN LISTS ARGN)
    if(NOT IS_ABSOLUTE "${s}")
      get_filename_component(s "${s}" ABSOLUTE)
    endif()
    list(APPEND SRCS "${s}")
  endforeach()

  # Append files in argument list to global property
  set_property(GLOBAL APPEND PROPERTY SRCS_LIST "${SRCS}")
endfunction(add_sources)