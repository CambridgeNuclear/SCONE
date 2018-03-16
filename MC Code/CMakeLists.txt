cmake_minimum_required(VERSION 2.8.9)
project(directory_test)

# Command to select Intel Fortran Compiler for testing purposes

#set(CMAKE_Fortran_COMPILER /opt/intel/Compiler/11.1/046/bin/intel64/ifort)

set(CMAKE_Fortran_COMPILER /home/mak60/gcc-5.2.0/bin/gfortran)

#set(CMAKE_CXX_COMPILE_ID "Intel")

enable_language(Fortran) 

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modFiles)
set(CMAKE_VERBOSE_MAKEFILE on)

# Add compilers flags
#if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
  set(CMAKE_Fortran_FLAGS  "-std=f2008 -pedantic -O3 -Wall") 
#endif()
# -flto = inlining at link time for now causes an error 

MESSAGE( STATUS CMAKE_fortran_COMPILER_ID)

##########################################################################
macro (add_sources)
    file (RELATIVE_PATH _relPath "${PROJECT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
    foreach (_src ${ARGN})
        if (_relPath)
            list (APPEND SRCS "${_relPath}/${_src}")
        else()
            list (APPEND SRCS "${_src}")
        endif()
    endforeach()
    if (_relPath)
        # propagate SRCS to parent directory
        set (SRCS ${SRCS} PARENT_SCOPE)
    endif()
endmacro()
########################################################################

SET (SRCS test.f90) 


add_library(DataStructures ./DataStructures/dictionary_class.f90
                           ./DataStructures/IOdictionary_class.f90)
target_link_libraries(DataStructures SharedModules ParticleObjects RandomNumbers)

# include all directories with libraries 
add_subdirectory(SharedModules)
add_subdirectory(RandomNumbers)
add_subdirectory(ParticleObjects)
add_subdirectory(NuclearData)
add_subdirectory(geometryObjects)
add_subdirectory(CollisionOperator)

### TEMP commands to complie the dictionary 
# Define Files to compile 


# Define dependencies of the library 

####

add_executable(test.out test.f90)
target_link_libraries(test.out SharedModules
                               RandomNumbers 
                               ParticleObjects 
                               NuclearDataStructures 
                               NuclearData 
                               emissionENDF 
                               angleLawENDF      
                               DataStructures 
                               CollisionOperator 
                               Geometry)
