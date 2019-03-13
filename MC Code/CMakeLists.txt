cmake_minimum_required(VERSION 3.0)
project(SCONE)
##############################################################################
# DEFINE COMPILATION OPTIONS 

option(BUILD_TESTS "If is set to ON will compile tests" ON)
option(LTO "Enables link-time optimisation" ON) 
option(COVERAGE "Collect Coverage Info" ON) 
option(DEBUG "Enable extra run-time checks" ON) 
 
set(CMAKE_VERBOSE_MAKEFILE off)

# Include local cmake modules. TODO: Test on WINDOWS!!! 
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake") 

###############################################################################
# CONFIGURE FORTRAN COMPILER 

enable_language(Fortran) 
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modFiles)

# Add compilers flags
#
# -Wno-surprising flag is enabled to suppress a buggy warning about final
# subroutine it should be periodicly removed to see if other surprising 
# warining are not raised This buggy warning should be repaired in GCC 7.0 
#
# -flto (Link Time Optimisation) is SUPER important for performance. 
#Allows to inline function across diffrent source files (and modules)  
#

# Set core Fortran Compiler flags -> Will be used in all files 
set(CMAKE_Fortran_FLAGS "-std=f2008 -O3 -g -pedantic -Werror -Wall -Wno-surprising -Wno-unused-dummy-argument")

# Add extra flags to global  
if(COVERAGE)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -coverage")
endif()
  
# Set extra flags for scone library 
if(LTO)
  list(APPEND scone_extra_flags "-flto" "-ffat-lto-objects")
endif() 

# Add Debugging flags to scone library
if(DEBUG)
  set(APPEND scone_extra_flags "-fcheck=bounds -fcheck=do -fcheck=mem -fcheck=pointer -Waliasing") 
endif()

###############################################################################
# CHECK FOR DEPENDENCIES 

# Add environmental variable to default search directories 
list(APPEND CMAKE_PREFIX_PATH $ENV{LAPACK_INSTALL})

find_package(LAPACK REQUIRED )


# Dependencies for BUILD_TESTS
if (BUILD_TESTS)
    # FIND PYTHON INTERPRETER 
    # Sets PYTHONINTERP_FOUND & PYTHON_EXECUTABLE
    find_package(PythonInterp REQUIRED)

    find_package(PFUNIT REQUIRED) 
    add_library(pFUnit STATIC IMPORTED)
    set_property(TARGET pFUnit PROPERTY IMPORTED_LOCATION ${PFUNIT_LIBRARIES}) 
      
endif()

###############################################################################
# IMPORT FUNCTIONS TO COLLECT SOURCE FILES 
include(add_sources) 
include(add_unit_tests)

###############################################################################
# COLLECT ALL SOURCE AND TEST FILES 

# Include Nested Directories  
add_subdirectory(RandomNumbers) 
add_subdirectory(LinearAlgebra)
add_subdirectory(SharedModules)    
add_subdirectory(VTK)  
add_subdirectory(ParticleObjects)  
add_subdirectory(NamedGrids)

add_subdirectory(NuclearData)      
add_subdirectory(GeometryObjects)  
add_subdirectory(Tallies)

add_subdirectory(CollisionOperator)
add_subdirectory(TransportOperator)

add_subdirectory(MOCObjects) 
add_subdirectory(UserInterface)

add_subdirectory(PhysicsPackages)
add_subdirectory(DataStructures)
            
###############################################################################
# Compile SCONE static library 

# Copy source files collected in add_sources to variable SRCS
get_property(SRCS GLOBAL PROPERTY SRCS_LIST)

# Compile library 
add_library(scone STATIC ${SRCS}) 
target_compile_options(scone PRIVATE ${scone_extra_flags} )
target_link_libraries(scone  ${LAPACK_LIBRARIES} )

###############################################################################
# COMPILE SOLVERS 
add_executable(test.out ./Apps/test_MAK.f90)
target_link_libraries(test.out scone )  
#add_executable(eigenCE.out        ./Apps/eigenCE.f90           ${SRCS}) 
#add_executable(eigenMG.out        ./Apps/eigenMG.f90           ${SRCS}) 
#add_executable(test_init_geom.out ./Apps/test_init_geom.f90    ${SRCS}) 
#add_executable(dictTest.out       ./Apps/dictTest.f90          ${SRCS})
#add_executable(parser.out         ./Apps/testACE.f90           ${SRCS})
#add_executable(eigenvalue.out     ./Apps/eigenvalue.f90        ${SRCS})
#add_executable(surfTests.f90       ./Apps/surfTests.f90      ${SRCS})
#add_executable(geomPrint.out       ./Apps/geomPrint.f90      ${SRCS})
#add_executable(rand.out        ./Apps/produceRNG.f90 ${SRCS})

add_executable(scone.out ./Apps/scone.f90 )
target_link_libraries(scone.out scone ) 

###############################################################################
# COMPILE TESTS 
if(BUILD_TESTS) 

  # Create directory in binary_dir for temporary preprocessed test files 
  file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/generated)
  
  # Create required test suites listing file  
  file(WRITE ${PROJECT_BINARY_DIR}/generated/testSuites.inc "")
 
  include_directories(
      ${PROJECT_BINARY_DIR}/generated
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
       COMMAND ${PYTHON_EXECUTABLE} ${PFUNIT_PREPROC}/pFUnitParser.py ${_testPath} ${PROJECT_BINARY_DIR}/generated/${testName}${extension}
       DEPENDS pFUnit ${_testPath}
       COMMENT "Preprocessing test ${testName}"
       VERBATIM
       )
    # APPEND LIST OF ALL PREPROCESSED UNIT TEST FILES & ADD TEST SUITE TO testSuites.inc 
    set(PREP_UNIT_TESTS ${PREP_UNIT_TESTS} ${PROJECT_BINARY_DIR}/generated/${testName}${extension})   
    file(APPEND ${PROJECT_BINARY_DIR}/generated/testSuites.inc "ADD_TEST_SUITE(${testName}_suite)\n")    
  endforeach()


  add_executable(allTests ${PFUNIT_INCLUDE_DIRS}/driver.F90  ${PREP_UNIT_TESTS}) 
  target_link_libraries(allTests pFUnit scone)  

  # Switch off warnings when compiling tests. Has some unused variables and invalid preprocessor 
  # directives  
  # Select GNU compiler and Build Robust Test Suite 
  set_target_properties(allTests PROPERTIES COMPILE_FLAGS "-w -DGNU -DBUILD_ROBUST")
  
endif() 