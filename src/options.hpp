#ifndef PT_OPTIONS_HPP_
#define PT_OPTIONS_HPP_

#include <string>

namespace pt {

struct Options {
  std::string newick_path;
  std::string fasta_path;
  std::string raxml_path;

  double lnl_offset;
  size_t thread_count;

  int optimization_radius;
  bool try_all_moves;

  bool skip_filtering;

  std::string good_trees_path;
  std::string visited_trees_path;
};

Options ParseArguments(int argc, const char* argv[]);

} // namespace pt


#endif // PT_OPTIONS_HPP_
