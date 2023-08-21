####################################################################################################
# CMAKE requirements and project INFO

cmake_minimum_required(VERSION 3.10)
project(SCONE)

# Set CMAKE Policies
#   CMP0074: Starting with CMake 3.12, all `find_package` commands will use `<package>_ROOT`
#            variables as hints for where to search. This policy provides compatibility
#            with older projects.
#            https://cmake.org/cmake/help/latest/policy/CMP0074.html
cmake_policy(SET CMP0074 NEW)

####################################################################################################
# DEFINE COMPILATION OPTIONS

option(BUILD_TESTS "If is set to ON will compile tests" ON)
option(LTO         "Enables link-time optimisation"     ON)
option(COVERAGE    "Collect Coverage Info"              OFF)
option(DEBUG       "Enable extra run-time checks"       OFF)
option(OPENMP      "Enable parallelism with OpenMP"     ON)

# Include local cmake modules. TODO: Test on WINDOWS!!!
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

####################################################################################################
# CONFIGURE FORTRAN COMPILER
#
# NOTE:
#   Compiler flags are set for gfortran
#

enable_language(Fortran)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modFiles)

# Add compilers flags
# The flags set below will be applied to all Fortran files
#
set(CMAKE_Fortran_FLAGS "-std=f2008 -O3 -g -pedantic -Wall -Wno-unused-dummy-argument -cpp")

# Add extra flags to all Fortran files
if(COVERAGE)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -coverage")
endif()

# Add Debugging flags to scone library
# NOTE: GFortran 7.5 may give buggy warnings here
#   If this is an issue, consider using `-fcheck=bounds,do,mem,pointer`
#   https://github.com/CambridgeNuclear/SCONE/pull/15#issuecomment-1384192425
if(DEBUG)
  list(APPEND scone_extra_flags "-fcheck=all,no-array-temps" "-Waliasing")
endif()


# If LTO is requested check if it is supported by compiler
# LTO flags will be set by CMAKE via `set_property` command later
if(LTO)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT has_lto OUTPUT lto_error)

  if(NOT has_lto)
      message(STATUS "<><><><><><><><><><><><><><><><><><><><><><><><><><>")
      message(FATAL_ERROR "LTO was requested but it is not supported by the compiler. "
                          "Use -DLTO=FALSE option to disable LTO. Error message: <${lto_error}>")
  endif()
endif()

# Set up OpenMP
if(OPENMP)
  find_package(OpenMP)

  if(NOT OpenMP_Fortran_FOUND)
    message(STATUS "<><><><><><><><><><><><><><><><><><><><><><><><><><>")
    message(FATAL_ERROR "It seems that the Fortran compiler does not support OpanMP"
                       "Use -DOPENMP=OFF to disable")
  endif()
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
endif()


####################################################################################################
# CHECK FOR DEPENDENCIES

# Add environmental variable to default search directories
list(APPEND CMAKE_PREFIX_PATH $ENV{LAPACK_INSTALL})

find_package(LAPACK REQUIRED )
message(STATUS ${LAPACK_LIBRARIES})

# Dependencies for BUILD_TESTS
if (BUILD_TESTS)
    # FIND PYTHON INTERPRETER
    # Sets PYTHONINTERP_FOUND & PYTHON_EXECUTABLE
    find_package(PythonInterp REQUIRED)

    find_package(PFUNIT 4.7 REQUIRED)
    add_library(pFUnit STATIC IMPORTED)
    set_property(TARGET pFUnit PROPERTY IMPORTED_LOCATION ${PFUNIT_LIBRARIES})

endif()

####################################################################################################
# IMPORT FUNCTIONS TO COLLECT SOURCE FILES
include(add_sources)
include(add_unit_tests)
include(add_integration_tests)

####################################################################################################
# COLLECT ALL SOURCE AND TEST FILES

# Include Nested Directories
add_subdirectory(RandomNumbers)
add_subdirectory(LinearAlgebra)
add_subdirectory(SharedModules)
add_subdirectory(Visualisation)
add_subdirectory(ParticleObjects)
add_subdirectory(NamedGrids)

add_subdirectory(NuclearData)
add_subdirectory(Geometry)
add_subdirectory(Tallies)

add_subdirectory(CollisionOperator)
add_subdirectory(TransportOperator)

add_subdirectory(UserInterface)

add_subdirectory(PhysicsPackages)
add_subdirectory(DataStructures)

####################################################################################################
# Compile SCONE static library

# Copy source files collected in add_sources to variable SRCS
get_property(SRCS GLOBAL PROPERTY SRCS_LIST)

# Compile library
add_library(scone STATIC ${SRCS})
target_compile_options(scone PRIVATE ${scone_extra_flags} )
target_link_libraries(scone  ${LAPACK_LIBRARIES} )

if(LTO)
  set_property(TARGET scone PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

####################################################################################################
# COMPILE SOLVERS
add_executable(scone.out ./Apps/scone.f90 )
target_link_libraries(scone.out scone )

####################################################################################################
# COMPILE UNIT TESTS
if(BUILD_TESTS)
  enable_testing()

  get_property(UNIT_TESTS GLOBAL PROPERTY UNIT_TESTS_LIST)

  # Make absolute paths relative to source directory
  foreach(_testPath IN LISTS UNIT_TESTS)
    file(RELATIVE_PATH _testPath ${CMAKE_SOURCE_DIR} ${_testPath})
    list(APPEND UNIT_TESTS_RELATIVE ${_testPath})
  endforeach()

  add_pfunit_ctest(unitTests
    TEST_SOURCES ${UNIT_TESTS_RELATIVE}
    LINK_LIBRARIES scone ${LAPACK_LIBRARIES}
    )

####################################################################################################
# COMPILE INTEGRATION TESTS

  get_property(INTEGRATION_TESTS GLOBAL PROPERTY INTEGRATION_TESTS_LIST)

  # Make absolute paths relative to source directory
  foreach(_testPath IN LISTS INTEGRATION_TESTS)
    file(RELATIVE_PATH _testPath ${CMAKE_SOURCE_DIR} ${_testPath})
    list(APPEND INTEGRATION_TESTS_RELATIVE ${_testPath})
  endforeach()

  add_pfunit_ctest(integrationTests
    TEST_SOURCES ${INTEGRATION_TESTS_RELATIVE}
    LINK_LIBRARIES scone ${LAPACK_LIBRARIES}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )

endif()

####################################################################################################
# Print summary of CMake options
#
message(STATUS "")
message(STATUS " ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*> ")
message(STATUS "  SCONE CONFIGURATION SUMMARY")
message(STATUS "  Unit Tests:                 " ${BUILD_TESTS})
message(STATUS "  Link-time optimisation:     " ${LTO})
message(STATUS "  Code coverage logging:      " ${COVERAGE})
message(STATUS "  Extra runtime debug checks: " ${DEBUG})
message(STATUS "  OpenMP parallelism:         " ${OPENMP})
message(STATUS "  Fortran compiler:     " ${CMAKE_Fortran_COMPILER})
message(STATUS "  Compiler version:     " ${CMAKE_Fortran_COMPILER_VERSION})
message(STATUS "  OpenMP version:       " ${OpenMP_CXX_VERSION})
message(STATUS "  Global Fortran flags: " ${CMAKE_Fortran_FLAGS})
message(STATUS "  SCONE library flags:  " "${scone_extra_flags}")
message(STATUS " <*))))><  <*))))><  <*))))><  <*))))><  <*))))><  <*))))>< ")
message(STATUS "")
