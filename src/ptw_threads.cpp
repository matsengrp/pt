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
#include "options.hpp"
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
  pt::Options options = pt::ParseArguments(argc, argv);

  //
  // load tree, sequences, and model info
  //

  pll_utree_t* tree = pll_utree_parse_newick(options.newick_path.c_str());

  std::vector<std::string> labels;
  std::vector<std::string> sequences;
  pt::pll::ParseFasta(options.fasta_path, tree->tip_count, labels, sequences);

  pt::pll::ModelParameters parameters = pt::pll::ParseRaxmlInfo(options.raxml_path);

  //
  // initialize the guru
  //

  pt::Guru guru(options.lnl_offset, options.thread_count, tree, parameters,
                labels, sequences, options.try_all_moves);

  //
  // go!
  //

  guru.Start();
  guru.Wait();

  //
  // write output
  //

  if (!options.skip_filtering) {
    // ensure that the good tree table is filtered based on the final threshold
    guru.FilterGoodTreeTable();
  }

  WriteTreeTable(guru.GetGoodTreeTable(), options.good_trees_path);

  if (!options.visited_trees_path.empty()) {
    WriteTreeTable(guru.GetVisitedTreeTable(), options.visited_trees_path);
  }

  //
  // clean up and return
  //

  pll_utree_destroy(tree, pt::pll::cb_erase_data);

  return 0;
}
