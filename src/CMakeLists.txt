set(PT_SRC
  authority.cpp
  model_parameters.cpp
  ordered-tree.cpp
  partition.cpp
  pll_partition.cpp
  pll-utils.cpp
  wanderer.cpp)

add_library(pt ${PT_SRC})
add_library(pt-static STATIC ${PT_SRC})

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

add_executable(ptw ptw.cpp)
target_link_libraries(ptw pt-static pll Threads::Threads)
