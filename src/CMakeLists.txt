set(PT_SRC
  authority.cpp
  compressed_tree.cpp
  guru.cpp
  label_dictionary.cpp
  move_tester.cpp
  move_tester/always.cpp
  move_tester/branch_neighborhood_optimizer.cpp
  move_tester/single_branch_optimizer.cpp
  options.cpp
  ordered_tree.cpp
  position.cpp
  wanderer.cpp)

add_library(pt_static STATIC ${PT_SRC})
add_dependencies(pt_static ptpll_static)
target_link_libraries(pt_static ptpll_static)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

add_executable(ptw_threads ptw_threads.cpp)
add_dependencies(ptw_threads pt_static ptpll_static)
target_link_libraries(ptw_threads pt_static ptpll_static Threads::Threads)
