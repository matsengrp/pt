set(PT_SRC
  authority.cpp
  guru.cpp
  ordered_tree.cpp
  wanderer.cpp)

add_library(pt_static STATIC ${PT_SRC})

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

add_executable(ptw ptw.cpp)
target_link_libraries(ptw pt_static ptpll_static pll Threads::Threads)

add_executable(ptw_threads ptw_threads.cpp)
target_link_libraries(ptw_threads pt_static ptpll_static pll Threads::Threads)
