project(pt CXX)
cmake_minimum_required(VERSION 3.6.2)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g -pg -Wall -Wextra -Wshadow -pedantic")

enable_testing()

include_directories(lib/libcuckoo/src)

add_subdirectory(src)
add_subdirectory(test)
