set(PT_SRC
  partition.cpp
  pll-utils.cpp)

add_library(pt ${PT_SRC})
add_library(pt-static STATIC ${PT_SRC})
