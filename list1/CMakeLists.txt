cmake_minimum_required( VERSION 3.17 )
project( list1 )


set(CMAKE_BUILD_TYPE "Release")

set( CMAKE_CXX_STANDARD 17 )
set( CMAKE_CXX_FLAGS_DEBUG "-g3 -Wall -Wextra -O0" )
set( CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=native -Wall -Wextra" ) #-DNDEBUG declares NDEBUG macro, ignoring assertions

include_directories("../utils")

find_package( OpenMP )
if (OPENMP_FOUND)
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp" )
    set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}" )
endif ()

add_executable(warm_up warm_up.cpp)
add_executable(percolation percolation.cpp)
