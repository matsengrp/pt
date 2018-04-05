#ifndef PT_OPTIONS_HPP_
#define PT_OPTIONS_HPP_

#include <memory>
#include <string>

namespace pt {

// forward declaration to resolve circular dependency
class MoveTester;

struct Options {
  std::string newick_path;
  std::string fasta_path;
  std::string raxml_path;

  double lnl_offset = 0.0;
  size_t thread_count = 1;
  unsigned int rate_categories = 4;
  unsigned int poll_ms = 100;
  int optimization_radius = -1;

  bool skip_filtering = false;
  bool optimize_models = false;
  bool map_mode = false;
  bool marginal_mode = false;

  std::shared_ptr<const MoveTester> move_tester;

  std::string good_trees_path;
  std::string visited_trees_path;
  std::string tested_trees_path;

  bool track_tested_trees = true;
};

Options ParseArguments(int argc, const char* argv[]);

} // namespace pt


#endif // PT_OPTIONS_HPP_
