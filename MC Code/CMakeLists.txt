cmake_minimum_required(VERSION 2.8.9)
project(directory_test)

enable_language(Fortran) 

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modFiles)

# include all directories with libraries 
add_subdirectory(SharedModules)
add_subdirectory(RandomNumbers)
add_subdirectory(ParticleObjects)
add_subdirectory(NuclearData)


# Add compilers flags
if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
  set(CMAKE_Fortran_FLAGS  "-std=f2003")
endif()


add_executable(test.out test.f90)
target_link_libraries(test.out SharedModules RandomNumbers ParticleObjects NuclearData)
