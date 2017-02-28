include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../src)

add_executable(test_pt test_pt.cpp)
target_link_libraries(test_pt pll pt-static)

add_executable(test_ptw test_ptw.cpp)
target_link_libraries(test_ptw pll pt-static)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test-data
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME test_pt COMMAND test_pt)
add_test(NAME test_ptw COMMAND test_ptw)
