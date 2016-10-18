set(PT_SRC
  partition.cpp
  pll-utils.cpp)

add_library(pt ${PT_SRC})
add_library(pt-static STATIC ${PT_SRC})

add_executable(test_pt test.cpp)
target_link_libraries(test_pt pll pt-static)
