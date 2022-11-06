cmake_minimum_required( VERSION 3.17 )
project( complex_systems )


set(CMAKE_BUILD_TYPE "Release")

set( CMAKE_CXX_STANDARD 20 )
set( CMAKE_CXX_FLAGS_DEBUG "-g3 -Wall -Wextra -O0" )
set( CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=native -Wall -Wextra" ) #-DNDEBUG declares NDEBUG macro, ignoring assertions

include_directories("utils")
include_directories("extern")

add_subdirectory(list1)
add_subdirectory(list2)