cmake_minimum_required(VERSION 3.0)
project(SCONE)
##############################################################################
# DEFINE COMPILATION OPTIONS 

option(BUILD_TESTS "If is set to ON will compile tests" ON)


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
# FLAGS FOR CODE COVERAGE        
# -coverage -Wno-surprising 
#
set(CMAKE_Fortran_FLAGS  "-std=f2008 -pedantic -g -O3 -Wall -Wno-surprising -Wno-unused-dummy-argument -flto")   

###############################################################################
# CHECK FOR DEPENDENCIES 

# Dependencies for BUILD_TESTS


if (BUILD_TESTS)
    # FIND PYTHON INTERPRETER 
    # Sets PYTHONINTERP_FOUND & PYTHON_EXECUTABLE
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
        add_library( pFUnit STATIC IMPORTED ) 
        set_property(TARGET pFUnit PROPERTY IMPORTED_LOCATION ${PFUNIT_DIR}/lib/libpfunit.a)        
    else()
        MESSAGE(FATAL_ERROR "PFUNIT_INSTALL DIRECTORY MUST BE SET")     
        include(ExternalProject)  

       # SET CMAKE ARGUMENTS FOR pFUnit Compilation with CMAKE
        set(pFUnit_CMAKE_ARGS
            -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
            -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external/pFUnit
            -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
            )

        ExternalProject_Add(pFUnit 
            PREFIX          ${PROJECT_BINARY_DIR}/external/pFUnit
            GIT_REPOSITORY  https://git.code.sf.net/p/pfunit/code 
            GIT_TAG         master
            CMAKE_ARGS      ${pFUnit_CMAKE_ARGS}
            )

        # Point to directory with pFUnit MOD files and store pFUnit root directory
        set(PFUNIT_DIR ${PROJECT_BINARY_DIR}/external/pFUnit)
    endif()
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
  endforeach()
  # append to global list
  set_property(GLOBAL APPEND PROPERTY SRCS_LIST "${SRCS}")  
endfunction(add_sources)

###############################################################################
# DEFINE FUNCTION TO COLLETC ALL UNIT TEST SOURCE FILES 
function(add_unit_tests) 
  get_property(is_defined GLOBAL PROPERTY UNIT_TEST_LIST DEFINED)
  if(NOT is_defined)
    define_property(GLOBAL PROPERTY UNIT_TEST_LIST
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
  endforeach()
  # Append to Global List 
  set_property(GLOBAL APPEND PROPERTY UNIT_TESTS_LIST "${TESTS}")    
  
endfunction(add_unit_tests) 

###############################################################################
# COLLECT ALL SOURCE AND TEST FILES 

add_sources( ./DataStructures/dictionary_class.f90
             ./DataStructures/IOdictionary_class.f90
             ./DataStructures/intMap_class.f90
             ./DataStructures/charMap_class.f90
             ./DataStructures/stack_class.f90
             ./DataStructures/dynArray_class.f90)
             
add_unit_tests( ./DataStructures/Tests/dictionary_test.f90
                ./DataStructures/Tests/intMap_test.f90
                ./DataStructures/Tests/charMap_test.f90
                ./DataStructures/Tests/dynArray_test.f90)
            
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
# COPY GLOBAL PROPERTY  SRC_LIST TO VARIABLE SRCS  
get_property(SRCS GLOBAL PROPERTY SRCS_LIST)


#foreach(s IN LISTS SRCS)
#  file(RELATIVE_PATH rs "${CMAKE_CURRENT_SOURCE_DIR}" "${s}")
#  string(REGEX REPLACE "r$" "" o "${CMAKE_CURRENT_BINARY_DIR}/${rs}")
#  add_custom_command(
#    OUTPUT "${o}"
#    COMMAND ${CMAKE_COMMAND} -E copy "${s}" "${o}"
#    DEPENDS "${s}"
#    COMMENT "Creating ${o}"
#    VERBATIM
#    )
#  list(APPEND PREP_SRCS ${o})
#endforeach()

###############################################################################
# COMPILE SOLVERS 
add_executable(test.out           ./Apps/test_MAK.f90          ${SRCS})
#add_executable(eigenCE.out        ./Apps/eigenCE.f90           ${SRCS}) 
#add_executable(eigenMG.out        ./Apps/eigenMG.f90           ${SRCS}) 
#add_executable(test_init_geom.out ./Apps/test_init_geom.f90    ${SRCS}) 
#add_executable(dictTest.out       ./Apps/dictTest.f90          ${SRCS})
#add_executable(parser.out         ./Apps/testACE.f90           ${SRCS})
#add_executable(eigenvalue.out     ./Apps/eigenvalue.f90        ${SRCS})
#add_executable(surfTests.f90       ./Apps/surfTests.f90      ${SRCS})
#add_executable(geomPrint.out       ./Apps/geomPrint.f90      ${SRCS})
add_executable(scone.out       ./Apps/scone.f90      ${SRCS})

###############################################################################
# COMPILE TESTS 
if(BUILD_TESTS) 

    # CREATE DIRECTORY IN BINARY_DIR FOR TEMPORARY FILES TO COMPILE TESTS  
    file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/generated)
    file(WRITE ${PROJECT_BINARY_DIR}/generated/testSuites.inc "")
 
    include_directories(
        ${PROJECT_BINARY_DIR}/generated
        ${PFUNIT_DIR}/mod
        )


    # PREPROCESS COLLECTED TEST FILES FROM GLOBAL PROPERTY INTO VARIABLE
    # COPY THE FILES TO THE BUILD FOLDER 
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
            COMMAND ${PYTHON_EXECUTABLE} ${PFUNIT_DIR}/bin/pFUnitParser.py ${_testPath} ${PROJECT_BINARY_DIR}/generated/${testName}${extension}
            DEPENDS pFUnit ${_testPath}
            COMMENT "Preprocessing test ${testName}"
            VERBATIM
            )
        # APPEND LIST OF ALL PREPROCESSED UNIT TEST FILES & ADD TEST SUITE TO testSuites.inc 
        set(PREP_UNIT_TESTS ${PREP_UNIT_TESTS} ${PROJECT_BINARY_DIR}/generated/${testName}${extension})   
        file(APPEND ${PROJECT_BINARY_DIR}/generated/testSuites.inc "ADD_TEST_SUITE(${testName}_suite)\n")    
    endforeach()


    add_executable(allTests ${PFUNIT_DIR}/include/driver.F90  ${PREP_UNIT_TESTS} ${SRCS}) 
    target_link_libraries(allTests pFUnit)  

    # Switch off warnings when compiling tests. Has some unused variables and invalid preprocessor 
    # directives  
    set_target_properties(allTests PROPERTIES COMPILE_FLAGS "-w -DGNU -DBUILD_ROBUST")
endif() 