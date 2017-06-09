#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <libpll/pll.h>
#include <model_parameters.hpp>
#include <pll_partition.hpp>
#include <pll_util.hpp>

#include "guru.hpp"
#include "wanderer.hpp"

void WriteTreeTable(pt::TreeTable& trees, const std::string& path)
{
  std::ofstream ofs(path);

  auto table = trees.lock_table();
  for (auto& item : table) {
    ofs << item.second << "\t" << item.first << "\n";
  }
}

int main(int argc, const char* argv[])
{
  if (argc != 8) {
    std::cerr << "usage:\n\n"
              << "ptw\n"
              << "\t<input_tree>\n"
              << "\t<input_alignment>\n"
              << "\t<input_raxml_info>\n"
              << "\t<lnl_offset>\n"
              << "\t<output_good_trees>\n"
              << "\t<output_visited_trees>\n"
              << "\t<thread_count>\n";
    return 1;
  }

  //
  // parse arguments
  //

  std::string newick_path(argv[1]);
  std::string fasta_path(argv[2]);
  std::string raxml_path(argv[3]);

  double lnl_offset = std::stod(argv[4]);

  std::string good_trees_path(argv[5]);
  std::string visited_trees_path(argv[6]);

  size_t thread_count = std::stoi(argv[7]);

  if (lnl_offset > 0.0) {
    std::cerr << "error: expected an lnl_offset < 0.\n";
    return 1;
  }

  if (thread_count < 1) {
    std::cerr << "error: expected a thread_count > 0.\n";
    return 1;
  }

  //
  // load tree, sequences, and model info
  //

  pll_utree_t* tree = pll_utree_parse_newick(newick_path.c_str());

  std::vector<std::string> labels;
  std::vector<std::string> sequences;
  pt::pll::ParseFasta(fasta_path, tree->tip_count, labels, sequences);

  pt::pll::ModelParameters parameters = pt::pll::ParseRaxmlInfo(raxml_path);

  //
  // initialize the guru
  //

  bool try_all_moves = true;
  pt::Guru guru(lnl_offset, thread_count, tree, parameters, labels, sequences,
                try_all_moves);

  //
  // go!
  //

  guru.Start();
  guru.Wait();

  //
  // write output
  //

  // ensure that the good tree table is filtered based on the final threshold
  guru.FilterGoodTreeTable();

  WriteTreeTable(guru.GetGoodTreeTable(), good_trees_path);
  WriteTreeTable(guru.GetVisitedTreeTable(), visited_trees_path);

  //
  // clean up and return
  //

  pll_utree_destroy(tree, pt::pll::cb_erase_data);

  return 0;
}
