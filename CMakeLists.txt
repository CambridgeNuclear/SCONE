####################################################################################################
# CMAKE requirements and project INFO

cmake_minimum_required(VERSION 3.10)
project(SCONE)

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
if(DEBUG)
  list(APPEND scone_extra_flags "-fcheck=bounds,do,mem,pointer" "-Waliasing")
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

    find_package(PFUNIT REQUIRED)
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

  # Create directory in binary_dir for temporary preprocessed test files
  file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/generated)

  # Create required test suites listing file
  file(WRITE ${PROJECT_BINARY_DIR}/generated/testSuites.inc "")

  include_directories(
      ${PFUNIT_MOD}
      )

  # Preprocess collected test files from global property into variable
  # Copy files to the build folder

  set(PREP_UNIT_TESTS)
  get_property(UNIT_TESTS GLOBAL PROPERTY UNIT_TESTS_LIST)

  foreach(_testPath IN LISTS UNIT_TESTS)
    # OBTAIN EXTENSION
    string(REGEX MATCH "\\.[a-zA-Z0-9]+$" extension ${_testPath})

    # REMOVE EXTENSION
    string(REGEX REPLACE "${extension}$" "" _testPath_temp ${_testPath})

    # OBTAIN FILE NAME WITHOUT EXTENSION
    string(REGEX MATCH "[a-zA-Z0-9_]+$" testName ${_testPath_temp})

    # ADD RULE PREPROCESS ALL TEST FILES TO A SINGLE FOLDER IN BINARY_DIR
    add_custom_command(
       OUTPUT ${PROJECT_BINARY_DIR}/generated/${testName}${extension}
       COMMAND ${PYTHON_EXECUTABLE} ${PFUNIT_PREPROC}/pFUnitParser.py ${_testPath}
                                    ${PROJECT_BINARY_DIR}/generated/${testName}${extension}
       DEPENDS pFUnit ${_testPath}
       COMMENT "Preprocessing test ${testName}"
       VERBATIM
       )
    # APPEND LIST OF ALL PREPROCESSED UNIT TEST FILES & ADD TEST SUITE TO testSuites.inc
    set(PREP_UNIT_TESTS ${PREP_UNIT_TESTS} ${PROJECT_BINARY_DIR}/generated/${testName}${extension})
    file(APPEND ${PROJECT_BINARY_DIR}/generated/testSuites.inc "ADD_TEST_SUITE(${testName}_suite)\n")
  endforeach()

  add_executable(unitTests ${PFUNIT_INCLUDE_DIRS}/driver.F90  ${PREP_UNIT_TESTS})
  target_link_libraries(unitTests pFUnit scone)
  target_include_directories(unitTests PUBLIC ${PROJECT_BINARY_DIR}/generated )

  # Switch off warnings when compiling tests. Has some unused variables and invalid preprocessor
  # directives
  # Select GNU compiler and Build Robust Test Suite
  set_target_properties(unitTests PROPERTIES COMPILE_FLAGS "-w -DGNU -DBUILD_ROBUST")

####################################################################################################
# COMPILE INTEGRATION TESTS

  # Create directory in binary_dir for temporary preprocessed test files
  file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/generatedInt)

  # Create required test suites listing file
  file(WRITE ${PROJECT_BINARY_DIR}/generatedInt/testSuites.inc "")

  # Preprocess collected test files from global property into variable
  # Copy files to the build folder

  set(PREP_INTEGRATION_TESTS)
  get_property(INTEGRATION_TESTS GLOBAL PROPERTY INTEGRATION_TESTS_LIST)

  foreach(_testPath IN LISTS INTEGRATION_TESTS)
    # OBTAIN EXTENSION
    string(REGEX MATCH "\\.[a-zA-Z0-9]+$" extension ${_testPath})

    # REMOVE EXTENSION
    string(REGEX REPLACE "${extension}$" "" _testPath_temp ${_testPath})

    # OBTAIN FILE NAME WITHOUT EXTENSION
    string(REGEX MATCH "[a-zA-Z0-9_]+$" testName ${_testPath_temp})

    # ADD RULE PREPROCESS ALL TEST FILES TO A SINGLE FOLDER IN BINARY_DIR
    add_custom_command(
       OUTPUT ${PROJECT_BINARY_DIR}/generatedInt/${testName}${extension}
       COMMAND ${PYTHON_EXECUTABLE} ${PFUNIT_PREPROC}/pFUnitParser.py ${_testPath}
                                    ${PROJECT_BINARY_DIR}/generatedInt/${testName}${extension}
       DEPENDS pFUnit ${_testPath}
       COMMENT "Preprocessing test ${testName}"
       VERBATIM
       )
    # APPEND LIST OF ALL PREPROCESSED UNIT TEST FILES & ADD TEST SUITE TO testSuites.inc
    set(PREP_INTEGRATION_TESTS ${PREP_INTEGRATION_TESTS} ${PROJECT_BINARY_DIR}/generatedInt/${testName}${extension})
    file(APPEND ${PROJECT_BINARY_DIR}/generatedInt/testSuites.inc "ADD_TEST_SUITE(${testName}_suite)\n")
  endforeach()


  add_executable(integrationTests ${PFUNIT_INCLUDE_DIRS}/driver.F90  ${PREP_INTEGRATION_TESTS})
  target_link_libraries(integrationTests pFUnit scone)
  target_include_directories(integrationTests PUBLIC ${PROJECT_BINARY_DIR}/generatedInt )


  # Switch off warnings when compiling tests. Has some unused variables and invalid preprocessor
  # directives
  # Select GNU compiler and Build Robust Test Suite
  set_target_properties(integrationTests PROPERTIES COMPILE_FLAGS "-w -DGNU -DBUILD_ROBUST")

  ##################################################################################################
  # Register the test suites with CMake
  # NOTE: Integration tests must be run from source directory
  enable_testing()
  add_test(NAME Unit_Tests COMMAND unitTests)
  add_test(NAME Integration_Tests COMMAND integrationTests WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
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
