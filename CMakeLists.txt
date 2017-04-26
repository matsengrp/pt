project(pt CXX)
cmake_minimum_required(VERSION 3.0)

# set required standard to C++11 for cmake versions prior to 3.1
# http://stackoverflow.com/questions/10851247/how-to-activate-c-11-in-cmake/25836953#25836953
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
  message(FATAL_ERROR "${CMAKE_CXX_COMPILER} does not support C++11.")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g -pg -Wall -Wextra -Wshadow -pedantic")

enable_testing()

include_directories(lib/libcuckoo/src)

add_subdirectory(src)
add_subdirectory(test)
