set(PT_SRC
  ordered-tree.cpp
  partition.cpp
  pll_partition.cpp
  pll-utils.cpp
  wanderer.cpp)

add_library(pt ${PT_SRC})
add_library(pt-static STATIC ${PT_SRC})
