set(PT_SRC
  authority.cpp
  ordered_tree.cpp
  wanderer.cpp)

add_library(pt_static STATIC ${PT_SRC})

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

add_executable(ptw ptw.cpp)
target_link_libraries(ptw pt_static ptpll_static pll Threads::Threads)
