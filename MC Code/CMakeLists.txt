cmake_minimum_required(VERSION 3.0)
project(SCONE)


###############################################################################
# CONFIGURE FORTRAN COMPILER 

enable_language(Fortran) 

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modFiles)
set(CMAKE_VERBOSE_MAKEFILE on)

# Add compilers flags
#
# -Wno-surprising flag is enabled to suppress a buggy warning about final subroutine 
# It should be periodicly removed to see if other surprising warining are not raised 
# This buggy warning should be repaired in GCC 7.0 
#
# -flto (Link Time Optimisation) is SUPER important for performance. Allows to inline function 
#  across diffrent source files (and modules)          
#
#
set(CMAKE_Fortran_FLAGS  "-std=f2008 -pedantic -g -O3 -Wall -Wno-surprising -Wno-unused-dummy-argument -flto")   

###############################################################################
# CHECK FOR DEPENDENCIES 

# FIND PYTHON INTERPRETER 
find_package(PythonInterp)

if(NOT PYTHONINTERP_FOUND) 
  MESSAGE(FATAL_ERROR "pFUnit Requires Python to be present") 
endif()    

# BUILD AND INSTALL pFUnit
# Check if PFUNIT_INSTALL enviromental variable exists. If it does it means user uses her/his 
# installation of pFUnit. If enviromental variable does not exists download and install pFUnit from 
# SourceForge.
#
if(DEFINED ENV{PFUNIT_INSTALL}) 
   set(PFUNIT_DIR $ENV{PFUNIT_INSTALL})

else()
    include(ExternalProject)  

   # SET CMAKE ARGUMENTS FOR pFUnit Compilation with CMAKE
    set(pFUnit_CMAKE_ARGS
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external/pFUnit
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
        )

    ExternalProject_Add(pfunit 
         PREFIX          ${PROJECT_BINARY_DIR}/external/pFUnit
         GIT_REPOSITORY  https://git.code.sf.net/p/pfunit/code 
         GIT_TAG         master
         CMAKE_ARGS      ${pFUnit_CMAKE_ARGS}
         )

    # Point to directory with pFUnit MOD files and store pFUnit root directory
    set(PFUNIT_DIR ${PROJECT_BINARY_DIR}/external/pFUnit)
endif()
###############################################################################
# DEFINE FUNCTION TO COLLECT ALL SOURCE FILES
function(add_sources)
  get_property(is_defined GLOBAL PROPERTY SRCS_LIST DEFINED)
  if(NOT is_defined)
    define_property(GLOBAL PROPERTY SRCS_LIST
      BRIEF_DOCS "List of source files"
      FULL_DOCS "List of source files to be compiled in one library")
  endif()
  # make absolute paths
  set(SRCS)
  foreach(s IN LISTS ARGN)
    if(NOT IS_ABSOLUTE "${s}")
      get_filename_component(s "${s}" ABSOLUTE)
    endif()
    list(APPEND SRCS "${s}")
    MESSAGE(STATUS ${SRCS})
  endforeach()
  # append to global list
  set_property(GLOBAL APPEND PROPERTY SRCS_LIST "${SRCS}")  
endfunction(add_sources)

###############################################################################
# DEFINE FUNCTION TO COLLETC ALL TEST SOURCE FILES 
function(add_tests) 
  get_property(is_defined GLOBAL PROPERTY TEST_LIST DEFINED)
  if(NOT is_defined)
    define_property(GLOBAL PROPERTY TEST_LIST
      BRIEF_DOCS "List of all test files" 
      FULL_DOCS "List of all pFUnit test suite files for preprocessing & compilation") 
  endif()
  # Make Absolute Paths 
  set(TESTS)
  foreach(t IN LISTS ARGN)
    if(NOT IS_ABSOLUTE "${t}")
      get_filename_component(t "${t}" ABSOLUTE)
    endif()
    list(APPEND TESTS "${t}")
    MESSAGE(STATUS ${TESTS})
  endforeach()
  # Append to Global List 
  set_property(GLOBAL APPEND PROPERTY TESTS_LIST "${TESTS}")    
  
endfunction(add_tests) 

###############################################################################
# COLLECT ALL SOURCE AND TEST FILES 

add_sources( ./DataStructures/dictionary_class.f90
             ./DataStructures/IOdictionary_class.f90
             ./DataStructures/maps_class.f90
             ./DataStructures/stack_class.f90)
add_tests( ./DataStructures/RANDOMFILE)
# Include Nested Directories  
add_subdirectory(SharedModules)    
add_subdirectory(RandomNumbers)  
add_subdirectory(VTK)  
add_subdirectory(ParticleObjects)  

add_subdirectory(NuclearData)      
add_subdirectory(GeometryObjects)  
add_subdirectory(Tallies)

add_subdirectory(CollisionOperator)
add_subdirectory(TransportOperator)

add_subdirectory(MOCObjects) 
add_subdirectory(UserInterface)

add_subdirectory(PhysicsPackages)

###############################################################################
# PREPROCESS COLLECTED SOURCES FROM GLOBAL PROPERTY INTO VARIABLE
# COPY THE FILES TO THE BUILD FOLDER 
set(PREP_SRCS)
get_property(SRCS GLOBAL PROPERTY SRCS_LIST)

foreach(s IN LISTS SRCS)
  file(RELATIVE_PATH rs "${CMAKE_CURRENT_SOURCE_DIR}" "${s}")
  string(REGEX REPLACE "r$" "" o "${CMAKE_CURRENT_BINARY_DIR}/${rs}")
  add_custom_command(
    OUTPUT "${o}"
    COMMAND ${CMAKE_COMMAND} -E copy "${s}" "${o}"
    DEPENDS "${s}"
    COMMENT "Creating ${o}"
    VERBATIM
    )
  list(APPEND PREP_SRCS ${o})
endforeach()

###############################################################################
# PREPROCESS COLLECTED TEST FILES FROM GLOBAL PROPERTY INTO VARIABLE
# COPY THE FILES TO THE BUILD FOLDER 
set(PREP_TESTS)
get_property(TESTS GLOBAL PROPERTY TESTS_LIST)

foreach(t IN LISTS TESTS)
  file(RELATIVE_PATH relPath "${CMAKE_CURRENT_SOURCE_DIR}" "${t}")
  string(REGEX REPLACE "r$" "" o "${CMAKE_CURRENT_BINARY_DIR}/${relPath}")
  add_custom_command(
    OUTPUT "${o}"
    COMMAND ${CMAKE_COMMAND} -E copy "${t}" "${o}"
    DEPENDS "${t}"
    COMMENT "Creating ${o}"
    VERBATIM
    )
  list(APPEND PREP_TESTS ${o})
endforeach()


###############################################################################
# COMPILE SOLVERS 
add_executable(test.out           ./Apps/test_MAK.f90          ${PREP_SRCS})
#add_executable(eigenCE.out        ./Apps/eigenCE.f90           ${PREP_SRCS}) 
#add_executable(eigenMG.out        ./Apps/eigenMG.f90           ${PREP_SRCS}) 
#add_executable(test_init_geom.out ./Apps/test_init_geom.f90    ${PREP_SRCS}) 
#add_executable(dictTest.out       ./Apps/dictTest.f90          ${PREP_SRCS})
#add_executable(parser.out         ./Apps/testACE.f90           ${PREP_SRCS})
add_executable(eigenvalue.out     ./Apps/eigenvalue.f90        ${PREP_SRCS})
#add_executable(surfTests.f90       ./Apps/surfTests.f90      ${PREP_SRCS})
add_executable(geomPrint.out       ./Apps/geomPrint.f90      ${PREP_SRCS})
add_executable(scone.out       ./Apps/scone.f90      ${PREP_SRCS})

###############################################################################
# COMPILE TESTS 
