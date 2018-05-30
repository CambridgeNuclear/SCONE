cmake_minimum_required(VERSION 2.8.11)
project(directory_test)


###############################################################################
# CONFIGURE ENVIROMENT 

#set(CMAKE_Fortran_COMPILER /home/mak60/gcc-6.3.0/bin/gfortran)
enable_language(Fortran) 

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modFiles)
set(CMAKE_VERBOSE_MAKEFILE on)

# Add compilers flags
#
# -Wno-surprising flag is enabled to suppress a buggy warning about final subroutine 
# It should be periodicly removed to see if other surprising warining are not raised 
# This buggy warning should be repaired in GCC 7.0 
#
set(CMAKE_Fortran_FLAGS  "-std=f2008 -pedantic -Ofast -Wall -Wno-surprising") # -flto  

# -flto = Link Time Optimisation - SUPER important for performance 
MESSAGE( STATUS CMAKE_fortran_COMPILER_ID)

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
# COLLECT ALL SOURCE FILES 

add_sources( ./DataStructures/dictionary_class.f90
             ./DataStructures/IOdictionary_class.f90)

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


###############################################################################
# PREPROCESS COLLECTED SOURCES FROM GLOBAL PROPERTY INTO VARIABLE
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
# COMPILE SOLVERS 
add_executable(test.out           ./Apps/test_MAK.f90          ${PREP_SRCS})
add_executable(eigenCE.out        ./Apps/eigenCE.f90           ${PREP_SRCS}) 
add_executable(eigenMG.out        ./Apps/eigenMG.f90           ${PREP_SRCS}) 
add_executable(test_init_geom.out ./Apps/test_init_geom.f90    ${PREP_SRCS}) 
add_executable(dictTest.out       ./Apps/dictTest.f90          ${PREP_SRCS})
add_executable(parser.out         ./Apps/testACE.f90           ${PREP_SRCS})
