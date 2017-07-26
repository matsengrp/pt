include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../src)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

add_executable(test_ptw test_ptw.cpp)
add_dependencies(test_ptw pt_static ptpll_static)
target_link_libraries(test_ptw pt_static ptpll_static pll Threads::Threads)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test-data
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME test_ptw COMMAND test_ptw)
