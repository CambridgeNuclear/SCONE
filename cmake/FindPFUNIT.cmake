# TODO: Sort out the licence

#[=======================================================================[.rst:

Find pFUnit
===========

Find the pFUnit Fortran Unit Test Framework

General guide to writing Find Modules in CMake is available here:
https://cmake.org/cmake/help/git-master/manual/cmake-developer.7.html#modules

Expects to find environmental variable PFUNIT_INSTALL that points to the installation
directory of pFUnit.

Result Variables
^^^^^^^^^^^^^^^^

Sets some of the usual findModule variables:

``PFUNIT_FOUND``
  True if PFUNIT was found.
``PFUNIT_LIBRARIES``
  Library to link against
``PFUNIT_INCLUDE_DIR``
  Directory with driver.F90 file

Sets some pFUnit Specific variables as well:
``PFUNIT_MOD``
  Location of the module directory of pFUnit
``PFUNIT_PREPROC``
  Location of pFUnit Python Preprocessor

#]=======================================================================]

# TODO: understand what these lines actually do
unset(PFUNIT_LIBRARY CACHE)
unset(PFUNIT_INCLUDE_DIR CACHE)
unset(PFUNIT_MODULES CACHE)
unset(PFUNIT_PREPROC CACHE)

# Find path to the pFUnit test driver programme
find_path(PFUNIT_INCLUDE_DIR
          NAMES driver.F90
          PATHS ENV PFUNIT_INCLUDE_DIR
                ENV PFUNIT_INSTALL
                ENV INCLUDE
          PATH_SUFFIXES include)

# Find path to the pFUnit library
find_library(PFUNIT_LIBRARY
          NAMES libpfunit.a libpfunit pfunit.a pfunit
          PATHS ENV PFUNIT_LIBRARY
                ENV PFUNIT_INSTALL
          PATH_SUFFIXES lib )

# Find path to the pFUnit module directory
find_path(PFUNIT_MODULES
          NAMES pfunit_mod pfunit_mod.mod
          PATHS ENV PFUNIT_MODULES
                ENV PFUNIT_INSTALL
          PATH_SUFFIXES mod)

# Find path to the pFnit Python preprocessor script
find_path(PFUNIT_PREPROC
          NAMES pFUnitParser.py
          PATHS ENV PATH
                ENV PFUNIT_INSTALL
          PATH_SUFFIXES bin)

# Support the REQUIRED and QUIET arguments and set PFUNIT_FOUND to True if found
set(failMSG "pFUnit unit test framework was not found. Set PFUNIT_INSTALL environmental \
variable to the root of the pFUnit installation directory. If pFUnit is not installed on your \
system, in README you will find instruction on how to obtain it and compile it")


include (FindPackageHandleStandardArgs)
find_package_handle_standard_args(PFUNIT ${failMSG} PFUNIT_LIBRARY
                                                    PFUNIT_INCLUDE_DIR
                                                    PFUNIT_MODULES
                                                    PFUNIT_PREPROC)

# Handle success and find failure
if(PFUNIT_FOUND)
  set(PFUNIT_LIBRARIES       ${PFUNIT_LIBRARY})
  set(PFUNIT_INCLUDE_DIRS    ${PFUNIT_INCLUDE_DIR})
  set(PFUNIT_MOD             ${PFUNIT_MODULES})

else()
  message(FATAL_ERROR "pFUnit unit test framework was not found. Set PFUNIT_INSTALL environmental
  variable to the root of the pFUnit installation directory. If pFUnit is not installed on your \
  system in README you will find instruction on how to obtain it and compile it" )

endif()

mark_as_advanced(PFUNIT_LIBRARY PFUNIT_INCLUDE_DIR PFUNIT_MODULES)


