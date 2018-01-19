#ifndef PT_OPTIONS_HPP_
#define PT_OPTIONS_HPP_

#include <memory>
#include <string>

#include "move_tester.hpp"

namespace pt {

struct Options {
  std::string newick_path;
  std::string fasta_path;
  std::string raxml_path;

  double lnl_offset = 0.0;
  size_t thread_count = 1;
  unsigned int rate_categories = 4;

  std::shared_ptr<const MoveTester> move_tester;

  bool skip_filtering = false;
  bool optimize_models = false;
  bool map_mode = false;

  std::string good_trees_path;
  std::string visited_trees_path;
  std::string tested_trees_path;

  bool track_tested_trees = true;
};

Options ParseArguments(int argc, const char* argv[]);

} // namespace pt


#endif // PT_OPTIONS_HPP_
